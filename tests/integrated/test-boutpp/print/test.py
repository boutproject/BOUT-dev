import boutpp as bc
import sys

bc.init("-d test")
if len(sys.argv) > 1:
    bc.print(*sys.argv[1:])
else:
    bc.print("We can print to the log from python 🎉")
