def read_limix(filepath):
    """
    Read the HDF5 limix file format.

    Parameters
    ----------
    filepath : str
        File path.

    Returns
    -------
    dict
        Phenotype and genotype.
    """
    Y = fetch(filepath, "phenotype/matrix")
    rows = _read_attrs(filepath, "phenotype/row_header")
    cols = _read_attrs(filepath, "phenotype/col_header")
    Y = _convert_headers(Y, rows, cols, "outcome")
    Y = Y.rename(sample_ID="sample")
    Y.name = "phenotype"

    G = fetch(filepath, "genotype/matrix")
    rows = _read_attrs(filepath, "genotype/row_header")
    cols = _read_attrs(filepath, "genotype/col_header")
    G = _convert_headers(G, rows, cols, "candidate")
    G = G.rename(sample_ID="sample")
    G.name = "genotype"

    return {"phenotype": Y, "genotype": G}


class fetcher(object):
    """
    Fetch datasets from HDF5 files.

    Parameters
    ----------
    filename : str
        Filename to an HDF5 file.

    Examples
    --------
    .. doctest::

        >>> from limix.io import hdf5
        >>> from limix import file_example
        >>> from limix.sh import extract
        >>>
        >>> with file_example("data.h5.bz2") as filepath:
        ...     filepath = extract(filepath, verbose=False)
        ...     with hdf5.fetcher(filepath) as df:
        ...         X = df.fetch('/group/dataset')
        ...         print('%.4f' % X[0, 0].compute())
        -0.0453
    """

    def __init__(self, filename):
        self._filename = filename

    def __enter__(self):
        import h5py

        self._f = h5py.File(self._filename, "r")
        return self

    def fetch(self, data_path):
        r"""Fetch a HDF5 dataset.

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

    def __exit__(self, *args):
        self._f.close()


def fetch(fp, path):
    """
    Fetches an array from hdf5 file.

    Parameters
    ----------
    fp : str
        HDF5 file path.
    path : str
        Path inside the HDF5 file.

    Returns
    -------
    numpy.ndarray
        Array read from the HDF5 dataset.
    """
    import h5py

    with h5py.File(fp, "r") as f:
        return f[path][:]


def _see(f_or_filepath, root_name="/", show_chunks=False):
    """
    Shows a human-friendly tree representation of the contents of a hdf5 file.

    Parameters
    ----------
    file_or_filepath : str, file
        HDF5 file path or a handler to an open one.
    root_name : str, optional
        Group path to be the tree root. Defaults to ``"/"``.
    show_chunks : bool, optional
        ``True`` to show chunks; ``False`` otherwise. Defaults to ``False``.

    Returns
    -------
    str
        String representation of the tree.
    """
    import h5py

    if isinstance(f_or_filepath, str):
        with h5py.File(f_or_filepath, "r") as f:
            _tree(f, root_name, show_chunks)
    else:
        _tree(f_or_filepath, root_name, show_chunks)


def _findnth(haystack, needle, n):
    parts = haystack.split(needle, n + 1)
    if len(parts) <= n + 1:
        return -1
    return len(haystack) - len(parts[-1]) - len(needle)


def _visititems(root, func, level=0, prefix=""):
    if root.name != "/":
        name = root.name
        eman = name[::-1]
        i1 = _findnth(eman, "/", level)
        name = "/" + eman[:i1][::-1]
        func(prefix + name, root)
    if not hasattr(root, "keys"):
        return
    for k in root.keys():
        if root.file == root[k].file:
            _visititems(root[k], func, level + 1, prefix)
        else:
            _visititems(root[k], func, 0, prefix + root.name)


def _tree(f, root_name="/", show_chunks=False):
    import h5py

    _names = []

    def get_names(name, obj):
        if isinstance(obj, h5py.Dataset):
            dtype = str(obj.dtype)
            shape = str(obj.shape)
            if show_chunks:
                chunks = str(obj.chunks)
                _names.append("%s [%s, %s, %s]" % (name[1:], dtype, shape, chunks))
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
        ns = n.split("/")
        add_to_node(root, ns)

    def child_iter(node):
        from numpy import argsort, asarray

        keys = list(node.children.keys())
        indices = argsort(keys)
        indices = asarray(indices)
        vals = list(node.children.values())
        return list(asarray(vals)[indices])

    import asciitree

    print(asciitree.draw_tree(root, child_iter))


def _read_attrs(filepath, path):
    from pandas import DataFrame
    import h5py

    with h5py.File(filepath, "r") as f:

        h = dict()
        for attr in f[path].keys():
            h[attr] = f[path + "/" + attr][()]

            if h[attr].dtype.kind == "S":
                h[attr] = h[attr].astype("U")

        return DataFrame.from_dict(h)


def _convert_headers(X, rows, cols, colname):
    from xarray import DataArray

    coords = {k: ("sample", rows[k]) for k in rows.keys()}
    coords.update({k: (colname, cols[k]) for k in cols.keys()})
    return DataArray(X, dims=["sample", colname], coords=coords)
