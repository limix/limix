from __future__ import unicode_literals

import re
import sys
from ast import literal_eval
from os import chdir, getcwd
from os.path import abspath, dirname, join

from setuptools import find_packages, setup

try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser

PY2 = sys.version_info[0] == 2


def _unicode_airlock(v):
    if isinstance(v, bytes):
        v = v.decode()
    return v


if PY2:
    string_types = unicode, str
else:
    string_types = bytes, str


class setup_folder(object):
    def __init__(self):
        self._old_path = None

    def __enter__(self):
        src_path = dirname(abspath(sys.argv[0]))
        self._old_path = getcwd()
        chdir(src_path)
        sys.path.insert(0, src_path)

    def __exit__(self, *_):
        del sys.path[0]
        chdir(self._old_path)


def get_init_metadata(metadata, name):
    expr = re.compile(r"__%s__ *=[^\"]*\"([^\"]*)\"" % name)
    prjname = metadata['packages'][0]
    data = open(join(prjname, "__init__.py")).read()
    return re.search(expr, data).group(1)


def make_list(metadata):
    return metadata.strip().split('\n')


def if_set_list(metadata, name):
    if name in metadata:
        metadata[name] = make_list(metadata[name])


def set_long_description(metadata):
    df = metadata['description_file']
    metadata['long_description'] = open(df).read()
    del metadata['description_file']


def convert_types(metadata):
    bools = ['True', 'False']
    for k in metadata.keys():
        v = _unicode_airlock(metadata[k])
        if isinstance(metadata[k], string_types) and v in bools:
            metadata[k] = v == 'True'


def setup_package():
    with setup_folder():

        config = ConfigParser()
        config.read('setup.cfg')

        metadata = dict(config.items('metadata'))
        metadata['packages'] = find_packages()
        metadata['platforms'] = literal_eval(metadata['platforms'])

        metadata['version'] = get_init_metadata(metadata, 'version')
        metadata['author'] = get_init_metadata(metadata, 'author')
        metadata['author_email'] = get_init_metadata(metadata, 'author_email')
        metadata['name'] = get_init_metadata(metadata, 'name')

        with open('requirements.txt') as f:
            metadata['install_requires'] = f.read().splitlines()

        with open('test-requirements.txt') as f:
            metadata['tests_require'] = f.read().splitlines()

        if_set_list(metadata, 'classifiers')
        if_set_list(metadata, 'keywords')
        if_set_list(metadata, 'cffi_modules')

        if 'extras_require' in metadata:
            metadata['extras_require'] = literal_eval(
                metadata['extras_require'])

        if 'console_scripts' in metadata:
            metadata['entry_points'] = {
                'console_scripts': make_list(metadata['console_scripts'])
            }
            del metadata['console_scripts']

        set_long_description(metadata)
        convert_types(metadata)

        setup(**metadata)


if __name__ == '__main__':
    setup_package()
