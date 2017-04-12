#!/usr/bin/env bash

pip install pypandoc
python -c "from pypandoc import download_pandoc as dp; dp(targetfolder='~/bin/');"
