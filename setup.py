import os
import sys
from os.path import join

from Cython.Distutils import build_ext
from setuptools import find_packages, setup
from setuptools.extension import Extension

try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
except (OSError, IOError, ImportError):
    long_description = open('README.md').read()


def _check_gcc_cpp11(cc_name):
    import subprocess
    try:
        cmd = cc_name + ' -E -dM -std=c++11 -x c++ /dev/null > /dev/null'
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        return False
    return True


class build_ext_subclass(build_ext):
    def build_extensions(self):
        import platform
        from distutils import sysconfig
        if (hasattr(self.compiler, 'compiler') and
                len(self.compiler.compiler) > 0):
            cc_name = self.compiler.compiler[0]
            stdcpp = '-std=c++11'
            if 'gcc' in cc_name and not _check_gcc_cpp11(cc_name):
                stdcpp = '-std=c++0x'
            for e in self.extensions:
                e.extra_compile_args.append(stdcpp)
                e.extra_compile_args.append('-Wno-deprecated-declarations')
                e.extra_compile_args.append('-Wno-unused-local-typedefs')
                e.extra_compile_args.append('-Wno-sign-compare')
                e.extra_compile_args.append('-Wno-self-assign')
                e.extra_compile_args.append('-Wno-macro-redefined')
                e.extra_compile_args.append('-Wno-unused-const-variable')

            conf_vars = sysconfig.get_config_vars()
            if 'MACOSX_DEPLOYMENT_TARGET' in conf_vars and len(conf_vars[
                    'MACOSX_DEPLOYMENT_TARGET']) > 0:
                _v1, _v2 = conf_vars['MACOSX_DEPLOYMENT_TARGET'].split('.')
                if int(_v1) == 10 and int(_v2) < 9:
                    stdcpp = '--stdlib=libc++'
                    for e in self.extensions:
                        e.extra_compile_args.append(stdcpp)
        build_ext.build_extensions(self)


def mac_workaround(compatible):
    import platform
    from distutils import sysconfig

    conf_vars = sysconfig.get_config_vars()
    if compatible:
        conf_vars['MACOSX_DEPLOYMENT_TARGET'] = '10.9'
        return
    vers = platform.mac_ver()[0].split('.')
    if len(vers) == 3:
        conf_vars['MACOSX_DEPLOYMENT_TARGET'] =\
            vers[0] + '.' + vers[1]
    else:
        conf_vars['MACOSX_DEPLOYMENT_TARGET'] = platform.mac_ver()[0]


def extra_compile_args():
    if sys.platform.startswith('win'):
        return []
    return [
        '-Wno-comment', '-Wno-unused-const-variable',
        '-Wno-overloaded-virtual', '-Wno-uninitialized',
        '-Wno-delete-non-virtual-dtor', '-Wunused-variable'
    ]


def core_extension(reswig):
    import numpy as np

    def globr(root, pattern):
        import fnmatch

        matches = []
        for root, _, filenames in os.walk(root):
            for filename in fnmatch.filter(filenames, pattern):
                matches.append(os.path.join(root, filename))

        return matches

    def swig_opts():
        return [
            '-c++', '-outdir', join('limix', 'deprecated'), '-I' + join('src')
        ]

    def nlopt_files():
        src = open(join('External', 'nlopt_src.files')).readlines()
        src = [join('External', 'nlopt', s).strip() for s in src]
        hdr = globr(join('External', 'nlopt'), '*/*.h')
        return (src, hdr)

    (src, hdr) = nlopt_files()
    src.extend(globr(join('src', 'limix'), '*.cpp'))
    hdr.extend(globr(join('src', 'limix'), '*.h'))

    incl = ['src', 'External', join('External', 'nlopt')]
    incl = [join(i) for i in incl]
    folder = join('External', 'nlopt')
    incl += [join(folder, f) for f in os.listdir(folder)]
    incl = [i for i in incl if os.path.isdir(i)]
    incl.extend([np.get_include()])

    wrap_file = join('src', 'interfaces', 'python', 'limix_wrap.cpp')
    i_file = join('src', 'interfaces', 'python', 'limix.i')

    if os.path.exists(wrap_file):
        src.append(wrap_file)
    else:
        src.append(i_file)

    depends = src + hdr

    ext = Extension(
        'limix.deprecated._core',
        src,
        include_dirs=incl,
        extra_compile_args=extra_compile_args(),
        swig_opts=swig_opts(),
        depends=depends)

    return ext


def ensemble_extension():
    import numpy as np

    src = [join('cython', 'lmm_forest', 'SplittingCore.pyx')]
    incl = [join('External'), np.get_include()]
    depends = src
    ext = Extension(
        'limix.ensemble.SplittingCore',
        src,
        language='c++',
        include_dirs=incl,
        extra_compile_args=extra_compile_args(),
        depends=depends)
    from Cython.Build import cythonize
    return cythonize(ext)


def setup_package(reswig, compatible):
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    if sys.platform == 'darwin':
        mac_workaround(compatible)

    needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
    pytest_runner = ['pytest-runner>=2.9'] if needs_pytest else []

    install_requires = [
        "scikit-learn", "pandas", "scipy", "h5py", "matplotlib"
    ]
    setup_requires = [] + pytest_runner
    tests_require = ['pytest']

    metadata = dict(
        name="limix",
        version="2.0.0",
        description="A flexible and fast mixed model toolbox.",
        long_description=long_description,
        keywords='linear mixed models, GWAS, QTL, ' +
        'Variance component modelling',
        maintainer="Limix Developers",
        author="Danilo Horta, Christoph Lippert, Paolo Casale, Oliver Stegle",
        author_email="stegle@ebi.ac.uk",
        maintainer_email="stegle@ebi.ac.uk",
        packages=find_packages(),
        include_package_data=True,
        url='http://pmbio.github.io/limix/',
        install_requires=install_requires,
        setup_requires=setup_requires,
        tests_require=tests_require,
        zip_safe=False,
        license='Apache License 2.0',
        ext_modules=[core_extension(reswig)] + ensemble_extension(),
        cmdclass=dict(build_ext=build_ext_subclass),
        entry_points={
            'console_scripts': [
                'limix_runner=limix.scripts.limix_runner:entry_point',
                'mtSet_postprocess=limix.scripts.mtSet_postprocess:entry_point',
                'mtSet_preprocess=limix.scripts.mtSet_preprocess:entry_point',
                'mtSet_simPheno=limix.scripts.mtSet_simPheno:entry_point',
                'mtSet_analyze=limix.scripts.mtSet_analyze:entry_point',
                'limix_converter=limix.scripts.limix_converter:entry_point',
                'iSet_analyze=limix.scripts.iSet_analyze:entry_point',
                'iSet_postprocess=limix.scripts.iSet_postprocess:entry_point',
            ]
        },
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Natural Language :: English',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ], )

    # http://stackoverflow.com/a/29634231
    import distutils.sysconfig
    cfg_vars = distutils.sysconfig.get_config_vars()
    for key, value in list(cfg_vars.items()):
        if type(value) == str:
            cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)


if __name__ == '__main__':
    reswig = False
    if "--reswig" in sys.argv:
        reswig = True
        sys.argv.remove("--reswig")

    compatible = False
    if "--compatible" in sys.argv:
        compatible = True
        sys.argv.remove("--compatible")

    setup_package(reswig, compatible)

