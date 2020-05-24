#!/usr/bin/env bash

function _awk () {
    awk '{
    if (match($0,"^$")) {
        empty=1
        #print "empty", $0
        if (running)
            out=1
        else
            out=0
        
    } else {
        if (match($0,"\\[ *RUN *\\]")) {
            #print "run", $0
            running=$0;
            out=0
        } else if (match($0,"\\[ *OK *\\]")) {
            #print "ok", $0
            running=""
            out=0
        } else if (match($0,"\\[ *FAILED *\\]")) {
            #print "failed", $0
            if (running)
                print running
            out=1
        } else if (match($0,"\\[-*\\]")) {
            #print "over", empty,"x",$0
            if (empty)
                out=0
            else
                out=1
        } else {
            if (running){
                print running
                running=""
            }
            out=1
        }
        empty=0
    }
    if (out)
        print $0
}'
}

set -o pipefail
for cmd in $@
do
    $cmd | _awk || exit
done
