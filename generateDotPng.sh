#!/bin/bash
DIR=/home/bmagalha/Workspace/neurox/build-qt/
FILES=$DIR/*.dot
cd $DIR
rm *.png
for f in $FILES
do
  dot -Tpng $f -o $f.png
done

