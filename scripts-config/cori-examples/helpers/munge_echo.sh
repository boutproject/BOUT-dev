#!/bin/bash

 

function munge() {

  case ":$1:" in

    *:$2:*) echo "$1" ;;

    ::) echo "$2" ;;

    *) echo "$2:$1" ;;

  esac

}

 

munged="$1"

 

shift

 

while [ -n "$1" ]; do

  munged="$(munge "$munged" "$1")"

  shift

done

 

echo "$munged"
