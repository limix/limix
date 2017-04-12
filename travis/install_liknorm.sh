#!/usr/bin/env bash

set -e

pushd .
git clone https://github.com/glimix/liknorm
cd liknorm
mkdir build
cd build
cmake ..
make
make test
sudo make install
sudo ldconfig
popd
