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
