#!/bin/bash

python Voffset.py

wdName=${PWD##*/}
rm list
ls -1tr vPlots/wd* > list


# while IFS= read -r line
# do
pdftk $(cat list) cat output vPlots/${wdName}."pdf"
# done < list
