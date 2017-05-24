import asciitree
import h5py


class h5data_fetcher(object):
    r"""
    Fetch datasets from HDF5 files.

    Parameters
    ----------
    filename : str
        Filename to an HDF5 file.

    Examples
    --------
    .. doctest::

        >>> from limix.io import h5data_fetcher
        >>> from limix.io.examples import hdf5_file_example
        >>> with h5data_fetcher(hdf5_file_example()) as df:
        ...     X = df.fetch('/group/dataset')
        ...     print('%.4f' % X[0, 0].compute())
        -0.0453
    """

    def __init__(self, filename):
        self._filename = filename

    def __enter__(self):
        self._f = h5py.File(self._filename, 'r')
        return self

    def fetch(self, data_path):
        r"""
        Fetch a HDF5 dataset.

        Parameters
        ----------
        data_path : str
            Path to a dataset.

        Returns
        -------
        X : dask array
        """
        from dask.array import from_array
        data = self._f[data_path]
        if data.chunks is None:
            chunks = data.shape
        else:
            chunks = data.chunks
        return from_array(data, chunks=chunks)

    def __exit__(self, *exc):
        self._f.close()


def fetch(fp, path):
    """Fetches an array from hdf5 file.
    :param str fp: hdf5 file path.
    :param str path: path inside the hdf5 file.
    :returns: An :class:`numpy.ndarray` representation of the corresponding hdf5 dataset.
    """
    with h5py.File(fp, 'r') as f:
        return f[path][:]


def see(f_or_filepath, root_name='/', ret=False, show_chunks=False):
    """Shows a human-friendly tree representation of the contents of
    a hdf5 file.
    :param f_or_filepath: hdf5 file path or a reference to an open one.
    :param str root_name: group to be the root of the tree.
    :param bool ret: Whether to return a string or print it.
    :param bool show_chunks: show the chunks.
    :returns str: String representation if is `ret=True`.
    """
    if isinstance(f_or_filepath, str):
        with h5py.File(f_or_filepath, 'r') as f:
            return _tree(f, root_name, ret, show_chunks)
    else:
        return _tree(f_or_filepath, root_name, ret, show_chunks)


def _findnth(haystack, needle, n):
    parts = haystack.split(needle, n + 1)
    if len(parts) <= n + 1:
        return -1
    return len(haystack) - len(parts[-1]) - len(needle)


def _visititems(root, func, level=0, prefix=''):
    if root.name != '/':
        name = root.name
        eman = name[::-1]
        i1 = _findnth(eman, '/', level)
        name = '/' + eman[:i1][::-1]
        func(prefix + name, root)
    if not hasattr(root, 'keys'):
        return
    for k in root.keys():
        if root.file == root[k].file:
            _visititems(root[k], func, level + 1, prefix)
        else:
            _visititems(root[k], func, 0, prefix + root.name)


def _tree(f, root_name='/', ret=False, show_chunks=False):
    _names = []

    def get_names(name, obj):
        if isinstance(obj, h5py.Dataset):
            dtype = str(obj.dtype)
            shape = str(obj.shape)
            if show_chunks:
                chunks = str(obj.chunks)
                _names.append("%s [%s, %s, %s]" % (name[1:], dtype, shape,
                                                   chunks))
            else:
                _names.append("%s [%s, %s]" % (name[1:], dtype, shape))
        else:
            _names.append(name[1:])

    _visititems(f, get_names)

    class Node(object):
        def __init__(self, name, children):
            self.name = name
            self.children = children

        def __str__(self):
            return self.name

    root = Node(root_name, dict())

    def add_to_node(node, ns):
        if len(ns) == 0:
            return
        if ns[0] not in node.children:
            node.children[ns[0]] = Node(ns[0], dict())
        add_to_node(node.children[ns[0]], ns[1:])

    _names = sorted(_names)
    for n in _names:
        ns = n.split('/')
        add_to_node(root, ns)

    def child_iter(node):
        from numpy import argsort, asarray

        keys = list(node.children.keys())
        indices = argsort(keys)
        indices = asarray(indices)
        vals = list(node.children.values())
        return list(asarray(vals)[indices])

    msg = asciitree.draw_tree(root, child_iter)
    if ret:
        return msg
    print(msg)
