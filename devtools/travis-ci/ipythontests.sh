#!/bin/sh
# Run ipython notebook tests

cd examples/fluorescence-binding-assay
testfail=0
shopt -s nullglob
for fn in *.ipynb; do
    echo "Testing IPython notebook $fn"
    ipynbtest.py "$fn" || testfail=1
done
cd ../..
if [ $testfail -eq 1 ]
then
    exit 1
fi

