#!/bin/bash

BASE_DIR="src/third_party/"

# Trova tutte le sottodirectory e itera su di esse
for dir in src/third_party/*/ ; do
    [ -L "{dir%/}" ] && continue
    cd "$dir"
    # echo "$dir"
    git checkout master && git pull
    cd -
done
