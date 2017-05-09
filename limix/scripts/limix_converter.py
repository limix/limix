import os
import time


def entry_point():
    from limix_legacy.deprecated.io.conversion import LIMIX_converter

    infostring = "limix_conveter.py, Copyright(c) 2014, "
    tim = time.ctime(os.path.getmtime(__file__))
    infostring += "The LIMIX developers\nlast modified: %s" % tim
    print(infostring)

    runner = LIMIX_converter(infostring=infostring)
    (options, args) = runner.parse_args()
    runner.run()
