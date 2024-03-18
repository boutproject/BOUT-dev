# /usr/bin/env python3

import jinja2
import os
import random
import string


def randomword(length):
    letters = string.ascii_lowercase
    return "".join(random.choice(letters) for i in range(length))


base, fn = __file__.rsplit("/", 1)

env = jinja2.Environment(loader=jinja2.FileSystemLoader(base))

template = env.get_template("field_tracking.hxx.jinja")

fn_tmp = f"{base}/.{fn}.{randomword(6)}"

with open(fn_tmp, "w") as f:
    f.write(template.render(ops="*/+-"))

os.rename(fn_tmp, __file__[:-3])
