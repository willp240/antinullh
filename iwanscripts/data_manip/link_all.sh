#!/bin/bash
dest=$1
shift
for filename in "$@"; do
    echo $filename
    ln -s $filename $dest
done
