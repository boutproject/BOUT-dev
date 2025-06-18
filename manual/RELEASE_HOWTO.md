# Making a New Release of BOUT++

This is checklist of things to do (in order) when making a new
release. This applies equally to both major/minor releases and bugfix
releases

- [ ] Check there are no open issues/PRs for this milestone
    - Fix or bump to next version
- [ ] Make a branch named `vX.Y.Z-rc`
    - The GitHub repo is setup to protect branches named in this style
    - Major and minor release points (`X`/`Y`) should be off
      `next`. Bugfix releases (`Z`) should be off `master`
- [ ] Make pull request into `master`
    - Any new bugfixes should be PRs into the RC branch, where
      "bugfixes" can include:
        - silencing warnings
        - improving documentation
        - adding tests
- [ ] Run `make check-all`
    - Raise issues for any tests that fail
- Possibly run `clang-tidy`, `clang-check`, `coverity`, etc.
- [ ] Review pinned pip package versions for CI
- [ ] Review bundled libraries
    
Before merging PR:

- [ ] Update locale translation files
    - `make -C locale update-all`
    - Be aware that this *will* update the timestamps and *possibly*
      reorder file paths in the .po and .pot files
- [ ] Update [`CHANGELOG.md`][changelog]:
    - Run [bout-changelog-generator.py LAST_RELEASE NEXT_RELEASE][bin/bout-changelog-generator.py]
        - See the docs for how to get the token
- [ ] Get list of authors:
    - [ ] `git log --format='%aN' | sort | uniq`
    - [ ] Compare to list in [`CITATION.cff`][citation], add new authors
- [ ] Prep a new Zenodo release:
    - https://doi.org/10.5281/zenodo.1423212
    - "New Version"
    - "Reserve DOI" -> copy DOI
    - Add any new authors
    - Save draft
- [ ] Change DOI in [`CITATION.cff`][citation] to new DOI
- [ ] Change DOI in [`README.md`][README] to new DOI
- [ ] Change date-released in [`CITATION.cff`][citation]
- [ ] Check `abidiff` to see if `soname` needs bumping in `makefile`:
- [ ] Change version number, removing prerelease tag in:
    - [ ]  [`configure.ac`][configure]: `AC_INIT`
    - [ ]  [`CITATION.cff`][citation]: `version`
    - [ ]  [`manual/sphinx/conf.py`][sphinx_conf]: `version` and `release`
    - [ ]  [`manual/doxygen/Doxyfile_readthedocs`][Doxyfile_readthedocs]: `PROJECT_NUMBER`
    - [ ]  [`manual/doxygen/Doxyfile`][Doxyfile]: `PROJECT_NUMBER`
    - [ ]  [`CMakeLists.txt`][CMakeLists]: `_bout_previous_version`, `_bout_next_version`
- [ ] Update what version of PETSc and SUNDIALS we support (upper bound)


After PR is merged:

- [ ] Make tarball: `make dist` from build directory. Ensure you are on a tag and correct version is used for archive and folder within.
- [ ] Try to summarise the changes!
- [ ] Make [GitHub Release][gh_release], include change summary **NB:** tag should have
      leading `v`
- [ ] Make tarball: `cmake -S . -B build && make dist -C build`
- [ ] Check tarball:
    - [ ] Is the folder name correct?
    - [ ] grep for the version - is always the released version used, not a pre-release?
- [ ] Upload tarball to GitHub Release
- [ ] Upload tarball to Zenodo and publish new version
- [ ] Email BOUT++ User Group mailing list, include change summary
- [ ] Make news post on project website, include change summary
- [ ] Update downloads page
- [ ] PR `master` into `next`
- [ ] Bump version number and add prerelease tag in:
    - [ ]  [`configure.ac`][configure]: `AC_INIT`
    - [ ]  [`CITATION.cff`][citation]: `version`
    - [ ]  [`manual/sphinx/conf.py`][sphinx_conf]: `version` and `release`
    - [ ]  [`manual/doxygen/Doxyfile_readthedocs`][Doxyfile_readthedocs]: `PROJECT_NUMBER`
    - [ ]  [`manual/doxygen/Doxyfile`][Doxyfile]: `PROJECT_NUMBER`
    - [ ]  [`CMakeLists.txt`][CMakeLists]: `_bout_previous_version`, `_bout_next_version`
    - [ ]  [`tools/pylib/_boutpp_build/backend.py`][backend.py]: `_bout_previous_version`, `_bout_next_version`

[Doxyfile]: ../manual/doxygen/Doxyfile
[Doxyfile_readthedocs]: ../manual/doxygen/Doxyfile_readthedocs
[citation]: ../CITATION.cff
[configure]: ../configure.ac
[sphinx_conf]: ../manual/sphinx/conf.py
[README]: ../README.md
[changelog]: ../CHANGELOG.md
[CMakeLists]: ../CMakeLists.txt
[backend.py]: ../tools/pylib/_boutpp_build/backend.py
[gcg]: https://github.com/github-changelog-generator/github-changelog-generator
[gh_release]: https://github.com/boutproject/BOUT-dev/releases/new
