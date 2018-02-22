test-initial
============

Test initial conditions are created correctly for varying number of
processes. Check all the initialisation expression functions
available:

|Name                                 | Description |
|-------------------------------------|--------------|
|abs(x)                               | Absolute value `|x|` |
|asin(x), acos(x), atan(x), atan(y,x) | Inverse trigonometric functions |
|ballooning(x)                        | Ballooning transform |
|ballooning(x,n)                      | Ballooning transform, using `n` terms (default 3) |
|cos(x)                               | Cosine |
|cosh(x)                              | Hyperbolic cosine |
|exp(x)                               | Exponential |
|tanh(x)                              | Hyperbolic tangent |
|gauss(x)                             | Gaussian `exp(-x^2/2) / sqrt{2pi}` |
|gauss(x, w)                          | Gaussian `exp[-x^2/(2w^2)] / (w sqrt{2pi})` |
|H(x)                                 | Heaviside function: `1` if `x > 0` otherwise `0` |
|log(x)                               | Natural logarithm |
|max(x,y,...)                         | Maximum (variable arguments) |
|min(x,y,...)                         | Minimum (variable arguments) |
|mixmode(x)                           | A mixture of Fourier modes |
|mixmode(x, seed)                     | seed determines random phase (default 0.5) |
|power(x,y)                           | Exponent `x^y` |
|sin(x)                               | Sine |
|sinh(x)                              | Hyperbolic sine |
|sqrt(x)                              | Square root |
|tan(x)                               | Tangent |
|erf(x)                               | The error function |
|TanhHat(x, width, centre, steepness) | The hat function |
|                                     | `0.5(tanh[s (x-[c-w/2])] + tanh[s (x-[c+w/2])] )` |

Unfortunately, the `ballooning` transform is a bit tricky to implement
in python, so this is currently skipped.

Another limitation is that it currently can't evaluate variable
references in the input file, so e.g.

    [foo]
    n = 4
    function = cos(foo:n * y)

won't work. This probably won't be fixed until a full expression
parser is implemented in python!

How it works
------------

`runtest` reads `data/BOUT.inp` and looks for sections that contain
the `function` attribute, ignoring three special ones, `var_x`,
`var_y`, and `var_z`, which are there simply to output the correct
coordinate grids.  Then it adds the following snippet of C++ to
`test_initial.cxx` for each of those sections that it finds:

    Field3D <name>;
    create_and_dump(<name>, "<name>");

where `<name>` is the name of each section. This creates, initialises
and writes out a `Field3D` according to the `function` attribute in
the corresponding section.

After running the test, `runtest` reads in the variables and checks
the answer against a python implementation (which in most cases is
just the Numpy implementation).
