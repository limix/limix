def file_type(filepath):
    if filepath.endswith('hdf5') or filepath.endswith('h5'):
        return 'hdf5'
    return 'unknown'
