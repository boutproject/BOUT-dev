Comments about the autogeneration
---------------------------------

The source code is generated using bash, which creates e.g. the `.pyx` files.
They in turn are compiled to `.cpp` files using cython.
Which are then compiled to machine code, with e.g. gcc.

Some comments about BASH
------------------------

Various stiles are used to print.
Here documents:
```
cat <<EOF
bla
EOF
```
In this here document variables are expanded. It is possible to call
functions, using ether the `\`function\`` or the `$(function)`
syntax. Thus care needs to be taken when writing documentation, as
they may contain backtics. Quotes must not be escaped.

There is also a version of the here document, which does not expand
```
cat <<"EOF"
bla
EOF
```
This version does not expand backticks, and also variables will not be
expanded.

echo statments exist again in two versions:
The version using double quotes `"` expands variables, and inline
functions are supported. echo statements using single quotes `'` are
not expanded.

`echo -n` does not print a new line after the statement.


The function `makelist` can be used in the fields section, to create
field specific strings, matching the dimensions and number of arguments.
