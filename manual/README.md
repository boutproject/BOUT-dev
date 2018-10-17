# BOUT++ documentation

The most up to date documentation consists of reStructuredText in the "sphinx"
subdirectory. This is used to create the online manual at
https://bout-dev.readthedocs.io/en/latest/. You can switch between the
manual for the last stable release and the latest development build by
clicking the "Read the Docs" bar in the left-hand sidebar, and
selecting either "stable" or "latest", respectively. You can also
download the manual as a PDF or as html from this menu too.

To build the manual locally, you need at least "sphinx" and
"recommonmark", which you can install using pip (or pip3):

```bash
$ pip install --user sphinx
$ pip install --user recommonmark
```

These documents can be built into a PDF using "sphinx-build":

```bash
$ make
```

To use e.g. "sphinx-build-3" instead of "sphinx-build", run
```bash
$ make sphinx-build=sphinx-build-3
```

This should create a file "BOUT.pdf" in the "manual" directory.

To get a local html version, run

```bash
$ make html
```

This should create a file "index.html" in the "manual/html" directory.

### API documentation

The majority of the codebase is documented using
[doxygen](www.doxygen.org). To build the API documentation, run

```bash
$ make doxygen
```

This creates html and LaTeX documentation under `doxygen/bout/html`
and `doxygen/bout/latex`.

It is possible to build the API documentation into the main manual
using "breathe". Install breathe:

```bash
$ pip install --user breathe
```

and comment out the following line in `sphinx/conf.py`:

```python
# Disable breathe
has_breathe = False
```

You can then build the sphinx documentation as normal.
