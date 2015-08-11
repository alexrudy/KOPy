#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pkg_resources
import cProfile
import os
import pyprof2html

def main():
    """Do the profiling work."""
    from KOPy.targets import TargetList
    filename = pkg_resources.resource_filename('KOPy.tests', 'data/big_starlist.txt')
    
    profile = cProfile.Profile()
    profile.runcall(TargetList.from_starlist, filename)
    
    output = os.path.relpath(os.path.join(os.path.dirname(__file__), "profiles", "targetlist_parse.profile"))
    try:
        os.makedirs(os.path.dirname(output) + "/")
    except IOError:
        pass
    
    profile.dump_stats(output)
    
    p2h = pyprof2html.Converter(output)
    outhtml = os.path.join(os.path.dirname(output), 'targetlist_parse')
    p2h.printout('html', outhtml)
    print("Profiled '{0:s}' in '{1:s}'".format('TargetList.from_starlist', outhtml))
    if p2h.profiledata_count > 20:
        p2h = pyprof2html.Converter(output)
        p2h.printout(filetype='html',
                     output_directory=outhtml,
                     output_htmlfile='index-all.html',
                     functions_number=99999)

if __name__ == '__main__':
    main()
    
    

