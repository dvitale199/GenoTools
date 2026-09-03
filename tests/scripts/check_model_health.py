#!/usr/bin/env python
# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Did this saved ancestry model converge, or did its training collapse?

Until round 19 the classifier's learning rate was never passed to XGBoost, so
gblinear descended at XGBoost's own default of 0.5 on its nondeterministic
Hogwild updater. Whether a fit converged or diverged was decided by a thread
race, and a diverged fit predicts one label for every sample while reporting
an accuracy exactly equal to that label's prevalence. The defect is present in
1.x and 2.0 alike, so every model already in production is worth checking.

The verdict comes from three numbers on the fitted booster -- largest |coef|,
largest |intercept|, and how many distinct classes it can emit -- which is
what cleanly separated the two real models behind the diagnosis:

    healthy GP2 model:   |coef|    1.03   |intercept| 4.24     10 classes
    collapsed PPMI model:|coef| 1727      |intercept| 3.0e15    1 class

Reading a 1.x model
-------------------
1.x pickled ``best_estimator_`` as a bare ``sklearn.pipeline.Pipeline``, and
those files unpickle in **no** environment available here: UMAP embeds numba
``Dispatcher`` objects carrying an ``impl_kind`` field that neither numba
0.63.1 nor 0.67.0 accepts, and sklearn objects a 1.3.0 file under 1.8.0. So
this reads them with a stub unpickler -- ``find_class`` hands back a
permissive dynamically-created *class* (not a function; ``NEWOBJ`` requires a
type) for every module outside a numpy/builtins allowlist. The
hyperparameters and fitted coefficients are plain scalars and arrays in the
instance ``__dict__``, so walking the object graph recovers them without any
of the code that wrote them.

Do not call ``repr()`` on a real cross-version sklearn estimator while doing
this: its pretty-printer touches attributes that did not exist in 1.3.0.

The coefficients are not in the ``__dict__``, though. ``XGBClassifier.coef_``
is derived on demand from the booster, and what gets pickled is
``Booster.handle`` -- a UBJSON buffer holding ``{"Config": ..., "Model": ...}``.
``Booster.load_model`` refuses that buffer (it is the internal config+model
pair, not ``save_model`` output), so this reads it with a small UBJSON parser
and takes ``learner.gradient_booster.model.weights``. For gblinear that is a
flat ``(num_feature + 1) x num_class`` array whose last row is the per-class
bias, which is where the 1e15 lives.

Usage
-----
    python tests/scripts/check_model_health.py <path> [<path> ...]

    # a 2.0 model directory
    python tests/scripts/check_model_health.py ~/.genotools/ref/my_model

    # a 1.x pickle
    python tests/scripts/check_model_health.py old_ancestry_model.pkl

Exit status is 1 if any model checked looks collapsed, so this can gate a
batch audit.
"""

import argparse
import io
import pickle
import struct
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

#: Above this a fitted linear model has diverged rather than converged. Kept
#: in step with `genotools.ancestry.fit_validation.MAX_HEALTHY_COEFFICIENT`,
#: and duplicated rather than imported so this script also runs against an
#: environment holding an older genotools.
MAX_HEALTHY_COEFFICIENT = 1e3

#: Modules whose classes are unpickled for real. Everything else is stubbed:
#: the numbers wanted here are numpy arrays and plain scalars, and the code
#: that produced them is exactly what will not import.
REAL_MODULES = ("numpy", "builtins", "collections", "copyreg", "_codecs")


#: UBJSON integer type markers, mapped to big-endian struct formats.
_UBJ_INT = {
    "i": (">b", 1),
    "U": (">B", 1),
    "I": (">h", 2),
    "l": (">i", 4),
    "L": (">q", 8),
}


class _Ubjson:
    """Just enough UBJSON to read an XGBoost booster buffer.

    Written out rather than taken from a library because the only thing needed
    is one weights array, and adding a dependency to a diagnostic script that
    must run against old environments is the wrong trade.
    """

    def __init__(self, buf: bytes):
        self.b, self.i = buf, 0

    def _u8(self) -> int:
        value = self.b[self.i]
        self.i += 1
        return value

    def _take(self, n: int) -> bytes:
        value = self.b[self.i : self.i + n]
        self.i += n
        return value

    def _int(self, marker: str) -> int:
        fmt, size = _UBJ_INT[marker]
        return struct.unpack(fmt, self._take(size))[0]

    def _key(self) -> str:
        return self._take(self._int(chr(self._u8()))).decode("utf-8", "replace")

    def value(self, marker: Optional[str] = None) -> Any:
        m = chr(self._u8()) if marker is None else marker
        if m == "Z":
            return None
        if m == "T":
            return True
        if m == "F":
            return False
        if m in _UBJ_INT:
            return self._int(m)
        if m == "d":
            return struct.unpack(">f", self._take(4))[0]
        if m == "D":
            return struct.unpack(">d", self._take(8))[0]
        if m == "S":
            return self._take(self._int(chr(self._u8()))).decode("utf-8", "replace")
        if m == "C":
            return self._take(1).decode("utf-8", "replace")
        if m == "{":
            return self._object()
        if m == "[":
            return self._array()
        raise ValueError(f"unhandled UBJSON marker {m!r} at byte {self.i}")

    def _header(self) -> Tuple[Optional[str], Optional[int]]:
        """The optional `$type` / `#count` prefix on a container."""
        element_type = count = None
        if chr(self.b[self.i]) == "$":
            self.i += 1
            element_type = chr(self._u8())
        if chr(self.b[self.i]) == "#":
            self.i += 1
            count = self._int(chr(self._u8()))
        return element_type, count

    def _object(self) -> Dict[str, Any]:
        element_type, count = self._header()
        out: Dict[str, Any] = {}
        # The key is bound to a local before the value is read: in
        # `out[self._key()] = self.value()` Python evaluates the value first,
        # which reads the stream out of order and mis-frames everything after.
        if count is not None:
            for _ in range(count):
                key = self._key()
                out[key] = self.value(element_type)
            return out
        while chr(self.b[self.i]) != "}":
            key = self._key()
            out[key] = self.value()
        self.i += 1
        return out

    def _array(self) -> List[Any]:
        element_type, count = self._header()
        out: List[Any] = []
        if count is not None:
            for _ in range(count):
                out.append(self.value(element_type))
            return out
        while chr(self.b[self.i]) != "]":
            out.append(self.value())
        self.i += 1
        return out


def booster_weights(handle: Any) -> Optional[Dict[str, Any]]:
    """Pull a gblinear booster's coefficients and biases out of its buffer.

    Returns None for anything that is not a readable gblinear model -- a tree
    booster carries `trees` rather than `weights`, and absence of a number
    here is not evidence of health.
    """
    try:
        payload = _Ubjson(bytes(handle)).value()
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None

    learner = payload.get("Model", payload)
    if isinstance(learner, dict):
        learner = learner.get("learner", learner)
    if not isinstance(learner, dict):
        return None

    model = learner.get("gradient_booster", {})
    model = model.get("model", {}) if isinstance(model, dict) else {}
    weights = model.get("weights") if isinstance(model, dict) else None
    if not weights:
        return None

    param = learner.get("learner_model_param", {}) or {}
    try:
        n_class = max(1, int(param.get("num_class", 1)))
        n_feature = int(param["num_feature"])
    except (KeyError, TypeError, ValueError):
        return None

    flat = np.asarray(weights, dtype=float)
    if flat.size != (n_feature + 1) * n_class:
        return None

    # gblinear lays the weights out as (num_feature + 1) x num_class, biases
    # last. That last row is where a diverged fit's 1e15 lives.
    table = flat.reshape(n_feature + 1, n_class)
    return {
        "max_abs_coefficient": _max_abs(table[:n_feature]),
        "max_abs_intercept": _max_abs(table[n_feature]),
        "n_classes": n_class,
        "n_features": n_feature,
        "boosted_rounds": model.get("boosted_rounds"),
    }


class _Stub:
    """Placeholder for a class from a module that will not import.

    Accepts any construction and any state, so the pickle machinery can build
    the object graph without the original code. The instance ``__dict__`` is
    the payload.
    """

    def __init__(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs

    def __setstate__(self, state):
        if isinstance(state, dict):
            self.__dict__.update(state)
        else:
            self._state = state

    def __reduce__(self):  # pragma: no cover - never re-pickled
        return (_Stub, ())

    # A few pickled objects are built by calling the class, others by
    # `__reduce__` with an iterator to extend or a mapping to update.
    def append(self, item):
        self.__dict__.setdefault("_items", []).append(item)

    def extend(self, items):
        self.__dict__.setdefault("_items", []).extend(items)

    def __setitem__(self, key, value):
        self.__dict__.setdefault("_mapping", {})[key] = value


class _StubUnpickler(pickle.Unpickler):
    """Unpickler that stubs out every class it cannot safely import."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.stubbed: set = set()

    def find_class(self, module, name):
        root = module.split(".")[0]
        if root in REAL_MODULES:
            return super().find_class(module, name)
        self.stubbed.add(f"{module}.{name}")
        # A *class*, not a function: the NEWOBJ opcode calls __class__.__new__.
        return type(f"Stub_{root}_{name}", (_Stub,), {"__module__": module})


def load_permissively(path: Path) -> Tuple[Any, set]:
    """Unpickle `path`, stubbing anything that will not import.

    A real load is tried first, because a 2.0 model written by this codebase
    unpickles cleanly and its enums and dataclasses then behave normally --
    stubbing an enum and reaching for a member of it raises.
    """
    raw = path.read_bytes()
    try:
        return pickle.loads(raw), set()
    except Exception:
        unpickler = _StubUnpickler(io.BytesIO(raw))
        return unpickler.load(), unpickler.stubbed


def walk(root: Any, max_nodes: int = 1_000_000):
    """Yield every object reachable from `root`, once each."""
    seen = set()
    stack = [root]
    while stack and len(seen) < max_nodes:
        node = stack.pop()
        marker = id(node)
        if marker in seen:
            continue
        seen.add(marker)
        yield node
        children: List[Any] = []
        state = getattr(node, "__dict__", None)
        if isinstance(state, dict):
            children += list(state.values())
        if isinstance(node, dict):
            children += list(node.keys()) + list(node.values())
        elif isinstance(node, (list, tuple, set, frozenset)):
            children += list(node)
        for child in children:
            if isinstance(child, (str, bytes, int, float, bool, type(None))):
                continue
            stack.append(child)


def _max_abs(value: Any) -> Optional[float]:
    """Largest absolute value in something array-shaped, or None."""
    try:
        array = np.abs(np.asarray(value, dtype=float))
    except (TypeError, ValueError):
        return None
    if array.size == 0:
        return None
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return float("inf")
    return float(finite.max())


def measure(path: Path) -> Dict[str, Any]:
    """Everything worth knowing about one saved model's numerical health."""
    target = path / "pipeline.pkl" if path.is_dir() else path
    if not target.exists():
        return {"path": str(path), "error": f"no such file: {target}"}

    try:
        obj, stubbed = load_permissively(target)
    except Exception as error:
        return {"path": str(path), "error": f"{type(error).__name__}: {error}"}

    coefs: List[float] = []
    intercepts: List[float] = []
    n_classes: Optional[int] = None
    params: Dict[str, Any] = {}
    metrics: Dict[str, Any] = {}

    for node in walk(obj):
        # A model that unpickled for real exposes `coef_` as a property
        # computed from the booster, so it never appears in `__dict__`. Read
        # it guardedly: `getattr` on an arbitrary node in the graph can run
        # someone else's `__getattr__` (six's lazy module importer will try to
        # import tkinter).
        for name, sink in (("coef_", coefs), ("intercept_", intercepts)):
            try:
                found = _max_abs(getattr(node, name))
            except Exception:
                continue
            if found is not None:
                sink.append(found)
        if n_classes is None:
            try:
                n_classes = int(np.asarray(node.classes_).size)
            except Exception:
                pass

        state = getattr(node, "__dict__", None)
        if not isinstance(state, dict):
            continue
        for key, value in state.items():
            if key == "handle":
                # An XGBoost Booster. `coef_` is derived on demand and never
                # pickled, so the numbers have to come out of this buffer.
                booster = booster_weights(value)
                if booster:
                    if booster["max_abs_coefficient"] is not None:
                        coefs.append(booster["max_abs_coefficient"])
                    if booster["max_abs_intercept"] is not None:
                        intercepts.append(booster["max_abs_intercept"])
                    if n_classes is None:
                        n_classes = booster["n_classes"]
                    params.setdefault("n_features", booster["n_features"])
                    if booster["boosted_rounds"] is not None:
                        params.setdefault(
                            "boosted_rounds", booster["boosted_rounds"]
                        )
            elif key == "coef_":
                found = _max_abs(value)
                if found is not None:
                    coefs.append(found)
            elif key == "intercept_":
                found = _max_abs(value)
                if found is not None:
                    intercepts.append(found)
            elif key == "classes_" and n_classes is None:
                try:
                    n_classes = int(np.asarray(value).size)
                except (TypeError, ValueError):
                    pass
            elif key in ("best_params_", "best_params"):
                if isinstance(value, dict):
                    params = {str(k): v for k, v in value.items()}
            elif key in ("kwargs", "_kwargs") and isinstance(value, dict):
                for name in ("lambda", "alpha", "learning_rate", "eta"):
                    if name in value:
                        params.setdefault(f"xgb__{name}", value[name])
            elif key in ("learning_rate", "n_estimators", "n_jobs", "nthread"):
                params.setdefault(key, value)
            elif key in ("train_accuracy", "test_accuracy", "best_score_"):
                metrics.setdefault(key, value)

    max_coef = max(coefs) if coefs else None
    max_intercept = max(intercepts) if intercepts else None
    diverged = any(
        value is not None and value > MAX_HEALTHY_COEFFICIENT
        for value in (max_coef, max_intercept)
    )

    # Every pre-round-19 model pickled `learning_rate=None`, because the field
    # was never passed. That makes it a reliable marker of which side of the
    # fix a saved model was trained on -- useful when auditing a directory of
    # them, since a converged pre-fix model is still evidence about luck rather
    # than about the process.
    pre_fix = "learning_rate" in params and params["learning_rate"] is None

    return {
        "path": str(path),
        "pickle": str(target),
        "pre_fix": pre_fix,
        "max_abs_coefficient": max_coef,
        "max_abs_intercept": max_intercept,
        "n_classes": n_classes,
        "diverged": diverged,
        "measurable": max_coef is not None or max_intercept is not None,
        "best_params": params,
        "reported_metrics": metrics,
        "stubbed_modules": sorted({name.split(".")[0] for name in stubbed}),
    }


def verdict(result: Dict[str, Any]) -> str:
    """One word for the top of the report."""
    if result.get("error"):
        return "UNREADABLE"
    if not result["measurable"]:
        return "UNMEASURABLE"
    if result["diverged"]:
        return "COLLAPSED"
    return "CONVERGED"


def describe(result: Dict[str, Any]) -> str:
    """The full report for one model."""
    label = verdict(result)
    lines = [f"{label}  {result['path']}"]
    if result.get("error"):
        lines.append(f"    {result['error']}")
        return "\n".join(lines)

    def show(value):
        return "n/a" if value is None else f"{value:.6g}"

    lines.append(
        f"    max |coef| {show(result['max_abs_coefficient'])}   "
        f"max |intercept| {show(result['max_abs_intercept'])}   "
        f"classes {result['n_classes'] if result['n_classes'] else 'n/a'}"
    )
    if result["best_params"]:
        shown = ", ".join(
            f"{key}={value}" for key, value in sorted(result["best_params"].items())
        )
        lines.append(f"    params: {shown}")
    if result["reported_metrics"]:
        shown = ", ".join(
            f"{key}={value}"
            for key, value in sorted(result["reported_metrics"].items())
        )
        lines.append(f"    reported: {shown}")
    if label == "COLLAPSED":
        lines.append(
            "    This model diverged during training. It predicts one label "
            "for every sample regardless of its reported accuracy; retrain it."
        )
    if label == "UNMEASURABLE":
        lines.append(
            "    No coefficients found. Absence is not health -- either the "
            "booster is a tree booster, or the pickle did not walk cleanly."
        )
    if result.get("pre_fix"):
        lines.append(
            "    Trained before the round-19 fix: the booster pickled "
            "learning_rate=None, so it descended at XGBoost's default of 0.5 "
            "on the nondeterministic Hogwild updater. This one came out fine, "
            "but the run that produced it was a coin flip."
        )
    if result["stubbed_modules"]:
        lines.append(
            f"    read with stubs for: {', '.join(result['stubbed_modules'])}"
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="2.0 model directories, or 1.x .pkl files",
    )
    args = parser.parse_args()

    collapsed = 0
    for path in args.paths:
        result = measure(path)
        print(describe(result))
        print()
        if verdict(result) == "COLLAPSED":
            collapsed += 1

    if collapsed:
        print(f"{collapsed} of {len(args.paths)} models checked have collapsed.")
        return 1
    print(f"{len(args.paths)} models checked, none collapsed.")
    return 0



if __name__ == "__main__":
    sys.exit(main())
