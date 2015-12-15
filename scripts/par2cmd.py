#!/bin/env python
from optparse import OptionParser
import sys, numpy, subprocess
'''
 converts a ftools par file into a command that can be used to execute the relevant FTOOL
'''
usage       = "Usage: %prog  [options] file.par"
description = "python script"
parser      = OptionParser(usage=usage,description=description)
parser.add_option("-r","--run"  ,default=0,action='store_true')
(opts, args)= parser.parse_args()
pfile = sys.argv[1]
tool = pfile.replace(".par","")+" "
pars = numpy.loadtxt(pfile,delimiter=',',dtype=str,unpack=True)
args = " ".join(["%s=%s"%(k,v) for k,v in zip(pars[0],pars[3])])
cmd = tool+args
print cmd
if opts.run:
    pts = subprocess.Popen(cmd,shell=True)
    pts.communicate()
