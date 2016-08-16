#!/bin/bash

python Vfit.py

wdName=${PWD##*/}
rm list
ls -1tr fitPlots/wd* > list


# while IFS= read -r line
# do
pdftk $(cat list) cat output fitPlots/${wdName}"_fit.pdf"
# done < list
