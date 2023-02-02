#!/bin/sh
# requires boutpp
# requires not make

set -ex

testit() {
    msg="$@"
    python3 test.py $msg > out.log
    grep "$msg" out.log
    grep "$msg" test/BOUT.log.0
}

testit Can we print to the log from python? 🎉  '__does_it_still__work {} {:s}'
