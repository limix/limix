#! /usr/bin/env python
# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
# All rights reserved.
#
# LIMIX is provided under a 2-clause BSD license.
# See LICENSE for the complete license.

import subprocess, sys, os
import time
import h5py
import scipy as sp
#import limix
from optparse import OptionParser
import pandas as pd
import limix_legacy.deprecated.io.output_writer as ow


class LIMIX_runner(object):
    '''
    A class for running a LIMIX pipeline.
    '''
    __slots__=["options","args","data", "result","infostring"]

    def __init__(self,infostring=None):
        '''
        nothing to initialize
        '''
        self.options=None
        self.args=None
        self.data=None
        self.result={}
        self.infostring=infostring
        if self.infostring is not None:
            self.result["infostring"]=self.infostring
        pass

    def parse_args(self):
        usage = "usage: %prog [options]"
        parser = OptionParser(usage=usage)
        parser.add_option("-D","--data_script", action="store", dest='data_script', type=str, help='A python script for loading the data', default=r"example_scripts/example_data.py")
        parser.add_option("-E","--experiment_script", action="store", dest='experiment_script', type=str, help='A python script for running the experiment', default=r"example_scripts/gwas_example.py")
        parser.add_option("-O","--outpath", action="store", dest='outpath', type=str, help='The output path (or file, if ending is .h5, the ooutput will be written to HDF5 file)', default="example_out")
        #parser.add_option('--seed', metavar='seed', type=int, help='The random seed', default=123123)
        parser.add_option("-T","--timestamp",action="store_true", dest='timestamp', help="Append unique timestamp to output value", default=False)
        parser.add_option("-R","--delimiter", action="store", dest='delimiter',type=str, help="The delimiter between output values", default="\t")
        parser.add_option("-F","--float_format", action="store", dest='float_format',type=str, help="Formating string for floating point output values", default="%.6e")
        (self.options, self.args) = parser.parse_args()
        self.result["options"]=str(self.options)
        return (self.options,self.args)

    def load_data(self):
        """
        Run the job specified in data_script
        """
        options=self.options
        command = open(self.options.data_script).read()
        self.result["data_script"]=command
        t0=time.time()
        data=None       #fallback data
        exec(command)   #creates variable data
        t1=time.time()
        print(("Elapsed time for data reading is %.2f seconds" % (t1-t0)))
        self.data=data
        return self.data

    def run_experiment(self):
        """
        Run the job specified in experiment_script
        """
        data=self.data
        options=self.options
        result=self.result

        command = open(self.options.experiment_script).read()
        result["experiment_script"]=command
        t0=time.time()
        exec(command)   #creates variable result
        t1=time.time()
        print(("Elapsed time for running the experiment is %.2f seconds" % (t1-t0)))
        self.result=result
        return self.result

    def write_resultfiles(self):
        """
        Write the output to disk
        """
        t0=time.time()
        writer = ow.output_writer(output_dictionary=self.result)
        if len(self.options.outpath)>=3 and self.options.outpath[-3:]==".h5":
            writer.write_hdf5(filename=self.options.outpath,timestamp=self.options.timestamp)
        else:
            writer.write_txt(outdir=self.options.outpath,timestamp=self.options.timestamp,delimiter=self.options.delimiter,float_format=self.options.float_format)
        t1=time.time()
        print(("Elapsed time for writing the output files is %.2f seconds" % (t1-t0)))


def entry_point():
    infostring = "limix_runner.py, Copyright(c) 2014, The LIMIX developers\nlast modified: %s" % time.ctime(os.path.getmtime(__file__))
    print (infostring)


    runner = LIMIX_runner(infostring=infostring)
    (options,args) = runner.parse_args()
    data = runner.load_data()
    result = runner.run_experiment()
    runner.write_resultfiles()
