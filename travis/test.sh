#!/bin/bash

set -e

python setup.py sdist

pip install dist/`ls dist | grep -i -E '\.(gz)$' | head -1` -vvv;

pushd /
python -c "import sys; import limix; sys.exit(limix.test())"
popd
