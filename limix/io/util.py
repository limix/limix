from os.path import exists


def file_type(filepath):
    imexts = ['.png', '.bmp', '.jpg', 'jpeg']
    if filepath.endswith('.hdf5') or filepath.endswith('.h5'):
        return 'hdf5'
    if filepath.endswith('.csv'):
        return 'csv'
    if filepath.endswith('.grm.raw'):
        return 'grm.raw'
    if filepath.endswith('.npy'):
        return 'npy'
    if _is_bed(filepath):
        return 'bed'
    if any([filepath.endswith(ext) for ext in imexts]):
        return 'image'
    return 'unknown'


def _is_bed(filepath):
    return all([exists(filepath + ext) for ext in ['.bed', '.bim', '.fam']])
