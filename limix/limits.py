import os


def is_large_file(filepath):
    large = 1024 * 1024 * 100
    return os.path.getsize(filepath) >= large
