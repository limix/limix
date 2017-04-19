#!/bin/bash

set -e

python setup.py sdist

if [ ${LEGACY} = "true" ]; then
    pip install dist/`ls dist | grep -i -E \'\.(gz)$\' | head -1`[legacy] -vvv;
else
    pip install dist/`ls dist | grep -i -E \'\.(gz)$\' | head -1` -vvv;
fi

pushd /
python -c "import sys; import limix; sys.exit(limix.test())"
popd
