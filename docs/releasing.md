# Cutting a release

There is no release automation in this repository — no publish workflow, no
release job. Publishing is a manual sequence run by a maintainer with PyPI
rights. This file is that sequence.

The version has one source of truth: `__version__` in `genotools/__init__.py`.
`setup.py` parses it, so nothing else needs editing.

---

## 1. Preflight

Run from a clean checkout of `main` with the release commit merged.

```bash
git checkout main && git pull
git status --short                       # must be empty of tracked changes
.venv/bin/python -m pytest tests/unit -q         # ~60s
.venv/bin/python -m pytest tests/regression -q   # ~7min, needs PLINK
```

Both suites must be green. The parity tests inside `tests/regression` skip
without `.venv-stable`; for a release, build it and run them for real:

```bash
bash tests/scripts/setup_stable_venv.sh v1.3.6
.venv/bin/python -m pytest tests/regression/test_parity.py -q
```

Do **not** set `PYTHONPATH` when running these. The 1.3.6 console script in
`.venv-stable` will import the 2.x working tree and fail for reasons unrelated
to the release. `cwd` alone is correct for both sides.

Then confirm the paperwork matches the code:

- [ ] `__version__` in `genotools/__init__.py` is the version being released
- [ ] `CHANGELOG.md` has an entry for it, and the entry describes what actually
      changed since the last **published** version — not since the last commit
- [ ] `MIGRATION_2.0.md` version references are current (it names specific
      versions in several places, including the `pip install` block and the
      model-provenance section)
- [ ] `docs/cli_args.md` documents every flag added or changed
- [ ] Any behavior that moves scientific results is stated plainly in both
      `CHANGELOG.md` and `MIGRATION_2.0.md`, including the ones that are
      *not* improvements
- [ ] **A pretrained model the released code can actually load is available**,
      and `README.md` names it correctly. `genotools-download` serves models
      from `https://storage.googleapis.com/genotools_refs/models/`, and the
      ones published there (`nba_v1`, `nba_v2`, `neurochip_v1`) are 1.x
      pickles that 2.x rejects by design. Until a 2.x-format model is published
      there, the documented getting-started path ends in a load error and the
      docs must say so, directing users to train their own with
      `--ref-panel`/`--ref-labels`

That last point is the one that costs users. A release note that says results
changed without saying whether the new results are better is honest; one that
implies an improvement nobody measured is not.

---

## 2. Build and inspect the artifact

Full detail, including why `build/` must be removed first, is in
[TESTING.md §8](../TESTING.md#8-building-a-release-artifact).

```bash
rm -rf build dist *.egg-info        # not optional
python -m pip wheel --no-deps -w dist .
```

Inspect before uploading:

```bash
python -c "import zipfile,glob; z=zipfile.ZipFile(glob.glob('dist/*.whl')[0]); \
  print(sum(i.file_size for i in z.infolist())/1024/1024, 'MB'); \
  print([i.filename for i in z.infolist() if not i.filename.endswith('.py')])"
```

A healthy 2.1.0 wheel is ~0.62 MB. If it is ~2.57 MB, it is shipping two 1.x
ancestry model pickles that the code cannot load — the `build/` directory was
stale.

Install it into a throwaway environment and check the version it reports:

```bash
python -m venv /tmp/relcheck && /tmp/relcheck/bin/pip install -q dist/*.whl
/tmp/relcheck/bin/python -c "import genotools; print(genotools.__version__)"
/tmp/relcheck/bin/genotools --help > /dev/null && echo "entry point ok"
```

(There is no `genotools --version` flag; the import above is the check. Adding
one would be a reasonable small change.)

---

## 3. Tag

Tag the exact commit the artifact was built from.

```bash
git tag -a v2.1.0 -m "GenoTools 2.1.0"
git push origin v2.1.0
```

The `v` prefix matters: `tests/scripts/setup_stable_venv.sh` and the CI parity
job resolve the baseline by tag name (`v1.3.6`).

---

## 4. Publish

```bash
python -m pip install --upgrade twine
twine check dist/*
twine upload dist/*                 # PyPI credentials required
```

Then verify from outside:

```bash
pip download --no-deps -d /tmp/pypicheck the_real_genotools==2.1.0
```

---

## 5. GitHub release

Create a release against the tag and paste the `CHANGELOG.md` section for this
version as the body. Attach nothing — PyPI holds the artifact.

---

## 6. After publishing

- [ ] Confirm the PyPI badge in `README.md` shows the new version
- [ ] If this release changes ancestry behavior, tell the people running
      production models — a release note reaches users who go looking, not
      users who have a model already trained and running
- [ ] Open the next development cycle by leaving `__version__` alone; it is
      bumped as part of the next release, not immediately after this one

---

## Worth automating later

Everything above is manual, which is why the 2.0 development versions were
bumped in `__init__.py` and never published — nothing connected the version
string to an actual release. A `release.yml` workflow triggered on tag push,
building the wheel and publishing via PyPI trusted publishing, would remove
both the credential handling and the "did anyone actually ship this?" question.
Not done here because publishing rights and the trusted-publisher configuration
are account-level decisions, not repository ones.
