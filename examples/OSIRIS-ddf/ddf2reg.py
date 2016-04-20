#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a .ddf file to a .reg file for DS9
"""

import os
from KOPy.instruments.osiris.regions import main

if __name__ == '__main__':
    dirname = os.path.dirname(__file__)
    starlist = os.path.join(dirname, "starlist.txt")
    ddf = os.path.join(dirname, "example.ddf")
    main((starlist, ddf, "--ds9", "--imdir", dirname))
    
