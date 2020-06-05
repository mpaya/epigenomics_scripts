#!/bin/bash

source config.txt

for f in ${exp[@]}; do mkdir -p "$basedir"/$f/{ARCHIVE,analysis}; done

[ -d "$refdir"] || mkdir -p "$refdir"

cp -r ../scripts "$basedir"

