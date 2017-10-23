from tempfile import mkdtemp
import random
import shutil
import os


class TmpDir(object):
    def __enter__(self):
        try:
            self.name = mkdtemp()
        except (PermissionError, FileExistsError):
            i = random.randint(0, 999999)
            tmp_name = ".tmp-limix-{}".format(i)
            try:
                os.mkdir(tmp_name)
                self.name = tmp_name
            except (PermissionError, FileExistsError):
                home = os.path.expanduser('~')
                tmp_name = os.path.join(home, tmp_name)
                os.mkdir(tmp_name)
                self.name = tmp_name
        return self.name

    def __exit__(self, *_):
        shutil.rmtree(self.name)
