from os.path import dirname, realpath, join
import tempfile
import shutil

_filenames = ["data.csv", "data.h5.bz2", "example", "pheno.csv"]


class file_example(object):
    def __init__(self, filename):
        if filename not in _filenames:
            raise ValueError(
                "Unrecognized file name. Choose one of these: {}".format(_filenames)
            )
        self._dirpath = tempfile.mkdtemp()
        self._filename = filename

    def __enter__(self):
        import pkg_resources

        filepath = join(self._dirpath, self._filename)

        if __name__ == "__main__":
            shutil.copy(join(dirname(realpath(__file__)), self._filename), filepath)
        else:
            resource_path = "example/{}".format(self._filename)
            content = pkg_resources.resource_string(__name__, resource_path)

            with open(filepath, "wb") as f:
                f.write(content)

        return filepath

    def __exit__(self, *_):
        shutil.rmtree(self._dirpath)
