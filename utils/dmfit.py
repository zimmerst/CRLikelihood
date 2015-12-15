#!/usr/bin/env python
""" 
Author: Alex Drlica-Wagner
Tiny script for translating DMFit channel
number into string and vice versa.
"""

INT2CHANNEL = {          
  2  :  "mu+mu-"       ,
  3  :  "tau+tau-"     ,
  4  :  "bb-bar"       ,
  5  :  "tt-bar"       ,
  6  :  "gluons"       ,
  7  :  "W+W-"         ,
  8  :  "Z Z"          ,
  9  :  "cc-bar"       ,
  10 :  "cosmo bb-bar" ,
  11 : "cosmo gam-gam" ,
}


def int2channel(int):
    return INT2CHANNEL[int]

CHANNEL2INT = {
   "mu+mu-"       : 2  ,
   "mumu"         : 2  ,
   "musrc"        : 2  ,
   "tau+tau-"     : 3  ,
   "tautau"       : 3  ,
   "tausrc"       : 3  ,
   "bb-bar"       : 4  ,
   "bb"           : 4  ,
   "bbbar"        : 4  ,
   "bbsrc"        : 4  ,
   "tt-bar"       : 5  ,
   "tt"           : 5  ,
   "gluons"       : 6  ,
   "W+W-"         : 7  ,
   "WW"           : 7  ,
   "wwsrc"        : 7  ,
   "Z Z"          : 8  ,
   "ZZ"           : 8  ,
   "cc-bar"       : 9  ,
   "cc"           : 9  ,
   "cosmo bb-bar" : 10 ,
   "cosmo gam-gam": 11 ,
}

def channel2int(channel):
    return CHANNEL2INT[channel]

def print_channels():
    for i in range(len(int2channel)):
        print i+2, int2channel[i+2]
        
if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()

    if len(args) == 1:
        arg = args[0]
        try:
            try:
                arg = int(arg)
                print int2channel[arg]
            except ValueError:
                print channel2int[arg]
        except KeyError:
            print "Model %s cannot be found"
            print_channels()
    else:
        print_channels()

