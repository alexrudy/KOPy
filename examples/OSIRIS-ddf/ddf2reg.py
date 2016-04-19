#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a .ddf file to a .reg file for DS9
"""

import os
from KOPy.instruments.osiris.ddf import DataDefinitionFile
from KOPy.instruments.osiris.regions import create_region_from_DDF
from astropy.coordinates import ICRS, SkyCoord
import astropy.units as u

DDF_FILE = os.path.join(os.path.dirname(__file__), "example.ddf")

if __name__ == '__main__':
    
    
    ddf = DataDefinitionFile.from_file(DDF_FILE)
    target = SkyCoord("10 59 18.11","+24 32 34.37", unit=(u.hourangle, u.degree))
    guidestar = SkyCoord("10 59 19.84","+24 32 45.80", unit=(u.hourangle, u.degree))
    target_frame = ddf.dataset.dithers.frame.at_origin(target)
    regions = create_region_from_DDF(ddf, target, guidestar)
    for region in regions():
        print(region)
