import xarray as xr
import numpy as np
import pandas as pd
import h5py
import limix
import os


def test_qtl_xarr(datadir):
    import pdb

    datadir.add("data/xarr.hdf5.gz")
    limix.sh.extract(os.path.join(datadir.tmpdir, "xarr.hdf5.gz"))
    pdb.set_trace()
    print()
    # limix.sh.extract
    # limix.io.hdf5.fetch()
