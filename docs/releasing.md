# Cutting a release

There is no release automation in this repository — no publish workflow, no
release job. Publishing is a manual sequence run by a maintainer with PyPI
rights. This file is that sequence.

The version has one source of truth: `__version__` in `genotools/__init__.py`.
`setup.py` parses it, so nothing else needs editing.

---

## Order of operations

Three publishes have to happen in this order, because each one is visible to
users the moment it lands:

1. **The model archive to GCS** (§2a). The released code's default is
   `nba_gp2_r12`; if PyPI goes first, every `genotools-download` between the two
   fails.
2. **The git tag** (§3), so the tag exists for the GitHub release and for
   `setup_stable_venv.sh` to resolve later.
3. **PyPI** (§4), last — it is the only irreversible step.

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
- [ ] **Every model named in `download_refs.MODELS` is actually uploaded**, and
      its recorded md5 matches the object in the bucket. The names are keys in
      that dict; the objects live at
      `gs://genotools_refs/models/<name>.zip`, and each archive must contain a
      single top-level directory matching its name, because `unzip_file` does a
      plain `extractall` into `<destination>/models`. A name in the dict with no
      object behind it fails at download; an object whose md5 has drifted fails
      checksum validation and exits 1. Verify with:

      ```bash
      python -c "from genotools.download_refs import MODELS; print(MODELS)"
      gsutil ls -L gs://genotools_refs/models/
      ```

That last point is the one that costs users. A release note that says results
changed without saying whether the new results are better is honest; one that
implies an improvement nobody measured is not.

---

## 2. Build and inspect the artifact

PyPI expects **both** a wheel and an sdist — every 1.x release shipped both — so
build with `build`, not `pip wheel`, which produces only the wheel. Full detail
on why `build/` must be removed first is in
[TESTING.md §8](../TESTING.md#8-building-a-release-artifact).

```bash
rm -rf build dist *.egg-info        # not optional
python -m pip install --upgrade build
python -m build                     # -> dist/*.whl and dist/*.tar.gz
```

Inspect both before uploading:

```bash
python - <<'EOF'
import glob, tarfile, zipfile
z = zipfile.ZipFile(glob.glob("dist/*.whl")[0])
print("wheel uncompressed:", round(sum(i.file_size for i in z.infolist())/1024/1024, 2), "MB")
print("wheel non-.py:", [i.filename for i in z.infolist() if not i.filename.endswith(".py")])
t = tarfile.open(glob.glob("dist/*.tar.gz")[0])
print("sdist entries:", len(t.getnames()))
print("pickles:", [n for n in z.namelist() + t.getnames() if n.endswith(".pkl")] or "none")
EOF
```

A healthy 2.1.0 wheel is ~0.62 MB uncompressed (~208 KB on disk) and the sdist
has ~78 entries. **Neither may contain a `.pkl`.** A ~2.57 MB wheel is shipping
two 1.x ancestry model pickles the code cannot load — the `build/` directory was
stale. Note the project has no `MANIFEST.in`, so the sdist's contents come from
setuptools' defaults plus `package_data`.

Install it into a throwaway environment and check the version it reports:

```bash
python -m venv /tmp/relcheck && /tmp/relcheck/bin/pip install -q dist/*.whl
/tmp/relcheck/bin/python -c "import genotools; print(genotools.__version__)"
/tmp/relcheck/bin/genotools --help > /dev/null && echo "entry point ok"
```

(There is no `genotools --version` flag; the import above is the check. Adding
one would be a reasonable small change.)

---

## 2a. Upload the model archive

Any model named in `download_refs.MODELS` must exist in the bucket with the
recorded md5:

```bash
gsutil cp nba_gp2_r12.zip gs://genotools_refs/models/
gsutil ls -l gs://genotools_refs/models/
```

Each archive must contain a single top-level directory matching its name —
`unzip_file` does a plain `extractall` into `<destination>/models`. Verify the
round trip before releasing:

```bash
genotools-download --model nba_gp2_r12 --destination /tmp/dlcheck
python -c "
from genotools.ancestry import AncestryModel
m = AncestryModel.load('/tmp/dlcheck/models/nba_gp2_r12')
print(len(m.common_snps), 'SNPs', list(m.label_encoder.classes_))"
```

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

**A PyPI upload is final.** A version number can never be reused, even after
deleting the release — if `2.1.0` goes up wrong, the only fix is `2.1.1`. Do
every check above first.

Credentials are an **API token**, not a password. On
[pypi.org](https://pypi.org) → *Account settings* → *API tokens*, create one
scoped to the `the-real-genotools` project; it is shown once and starts with
`pypi-`. The username is the literal string `__token__`.

```bash
python -m pip install --upgrade twine
twine check dist/*                  # metadata renders; run before every upload
twine upload dist/*                 # username: __token__   password: pypi-...
```

To avoid pasting the token each time, put it in `~/.pypirc` (mode 600):

```ini
[pypi]
  username = __token__
  password = pypi-AgEIcHlwaS5vcmc...
```

**Rehearse on TestPyPI first** if you want the upload path exercised without
consequences. It is a separate site with its own account and token, and it burns
the version number there too:

```bash
twine upload --repository testpypi dist/*
pip install --index-url https://test.pypi.org/simple/ \
  --extra-index-url https://pypi.org/simple the_real_genotools==2.1.0
```

Then verify the real thing from outside the machine that built it:

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
