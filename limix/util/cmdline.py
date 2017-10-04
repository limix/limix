import subprocess


def run_commandline(cmd):
    return subprocess.check_output(cmd, shell=True).decode()
