from os.path import dirname, realpath, join


def hdf5_file_example():
    dir_path = dirname(realpath(__file__))
    return join(dir_path, 'data', 'hdf5', 'data.h5')


def csv_file_example():
    dir_path = dirname(realpath(__file__))
    return join(dir_path, 'data', 'csv', 'data.csv')


def gen_file_example():
    dir_path = dirname(realpath(__file__))
    return join(dir_path, 'data', 'gen', 'example')


def numpy_kinship_file_example():
    dir_path = dirname(realpath(__file__))
    return join(dir_path, 'data', 'numpy', '1000G_kinship.npy')
