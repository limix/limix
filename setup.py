from setuptools import setup

if __name__ == '__main__':
    readme = open('README.md').read()
    try:
        import pypandoc
        long_description = pypandoc.convert_text(
            readme, 'rst', format='markdown')
    except (ImportError, RuntimeError, OSError):
        long_description = readme

    console_scripts = ["limix = limix.cmdlimix:entry"]
    setup(
        entry_points=dict(console_scripts=console_scripts),
        long_description=long_description)
