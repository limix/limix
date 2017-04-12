from dask.dataframe import read_csv as _read_csv

def read_csv(filename):
    df = _read_csv(filename)
    df.set_index(df.columns[0], inplace=True)
    return df
