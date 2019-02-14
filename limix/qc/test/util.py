from numpy.testing import assert_allclose, assert_, assert_array_equal


def _assert_attr_values(X, R):
    for attr in ["index", "columns"]:
        if hasattr(X, attr):
            assert_array_equal(getattr(X, attr), getattr(R, attr))


def assert_mat_proc(func, data, astype):
    def get():
        return astype(data["X"], data["samples"], data["candidates"]).copy()

    X = get()
    assert_allclose(func(X), data["R"])
    assert_allclose(X, get())
    assert_(isinstance(func(X), type(get())))
    _assert_attr_values(X, get())

    X = get()
    assert_allclose(func(X, axis=0), data["Rt"])
    assert_allclose(X, get())
    assert_(isinstance(func(X, axis=0), type(get())))
    _assert_attr_values(X, get())


def assert_mat_proc_inplace(func, data, astype):
    def get():
        return astype(data["X"], data["samples"], data["candidates"]).copy()

    X = get().copy()
    assert_allclose(func(X, inplace=True), data["R"])
    assert_allclose(X, data["R"])
    assert_(isinstance(func(X, inplace=True), type(get())))
    _assert_attr_values(X, get())

    X = get().copy()
    assert_allclose(func(X, axis=0, inplace=True), data["Rt"])
    assert_allclose(X, data["Rt"])
    assert_(isinstance(func(X, axis=0, inplace=True), type(get())))
    _assert_attr_values(X, get())
