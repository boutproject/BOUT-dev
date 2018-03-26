# BOUT++ Library - Write fluid simulations in curviilinear geometry
# Copyright (C) 2018 David Schwörer
#
# Contact: Ben Dudson, bd512@york.ac.uk
#
# This file is part of BOUT++.
#
# BOUT++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BOUT++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.


from common import *

test = UniqueList()
test.append(3)
test.append(4)
test.append(5)

ASSERT(str(test) == "[3, 4, 5]")

try:
    test.append(4)
    raise BaseException()
except RuntimeError:
    pass
