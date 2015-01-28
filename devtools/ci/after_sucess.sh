# Coveralls is temporarily commented out until we can set this up.
#coveralls

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


if [[ "2.7 3.3 3.4" =~ "$python" ]]; then
    echo "Uploading to binstar..."
    conda install --yes binstar
    echo binstar upload -t $BINSTAR_TOKEN --force -u omnia -p assaytools-dev $HOME/miniconda/conda-bld/linux-64/assaytools-dev-*
    binstar upload -t $BINSTAR_TOKEN --force -u omnia -p assaytools-dev $HOME/miniconda/conda-bld/linux-64/assaytools-dev-*
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi


# Create the docs and push them to S3
# -----------------------------------

conda install --yes pip
conda config --add channels http://conda.binstar.org/omnia
conda install --yes `conda build devtools/conda-recipe --output`
pip install numpydoc s3cmd msmb_theme
conda install --yes `cat docs/requirements.txt | xargs`

conda list -e

(cd docs && make html)
#python devtools/ci/push-docs-to-s3.py
#python devtools/ci/update-versions.py
