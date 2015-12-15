#!/usr/bin/env python
import time

LEAP_SECS = [
    time.mktime(time.strptime('2005-12-31 23:59:59', "%Y-%m-%j %H:%M:%S")),
    time.mktime(time.strptime('2009-12-31 23:59:59', "%Y-%m-%j %H:%M:%S")),
    time.mktime(time.strptime('2012-06-30 23:59:59', "%Y-%m-%j %H:%M:%S")),
]

def utc2met(t = None):
    if t is None:
        unixSecs=time.mktime(time.gmtime())
    else:
        unixSecs=time.mktime(time.strptime(t,"%Y-%m-%d %H:%M:%S"))
    
    # Fermi Mission Elapsed Time starts at 00:00:00 UTC Jan 1, 2001
    missionEpoch = time.mktime(time.strptime('Mon Jan 1 00:00:00 2001'))
    met = unixSecs - missionEpoch

    # Leap seconds
    for sec in LEAP_SECS:
        if unixSecs > sec: met += 1
    
    return met

if __name__ == "__main__":
    import sys
    from optparse import OptionParser
    usage = "Usage: %prog  [YYYY-MM-DD HH:MM:SS]"
    description = "Convert UTC to Fermi-LAT MET"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()

    if len(sys.argv) == 1:
        met = utc2met()
    else:
        met = utc2met(' '.join(sys.argv[1:]))

    print "%.0f"%(met)
