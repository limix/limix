import sys
import tempfile

PY2 = sys.version_info < (3, )

if PY2:
    import shutil
    from contextlib import contextmanager

    @contextmanager
    class TemporaryDirectory(object):
        def __init__(self):
            self._filepath = None

        def __enter__(self):
            self._filepath = tempfile.mkdtemp()
            return self._filepath

        def __exit__(self, *args):
            if self._filepath is not None:
                shutil.rmtree(self._filepath)
else:
    from tempfile import TemporaryDirectory
