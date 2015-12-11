#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make a UKIRT FS Starlist
"""

import collections
import numpy as np
from KOPy.targets import Target, TargetList
from astropy.coordinates import SkyCoord
import astropy.units as u
import datetime

VIZIER_UKIRT_FS_CATALOG = "J/MNRAS/373/781"
VIZIER_ADS_REFERENCE = "2006MNRAS.373..781L"

def row_to_target(row):
    """Make a target object out of a row from the Vizier UKIRT FS catalog."""
    
    # Get the position from Vizier. Use the Vizier normalized positions, which are in degrees J2000 (FK5).
    position = SkyCoord(row['_RAJ2000'] * u.deg, row['_DEJ2000'] * u.deg, frame='fk5')
    
    # Collect additional target keywords which might be useful from the Vizier catalog table.
    keywords = collections.OrderedDict()
    for filter in "JHK":
        # Collect magnitudes and errors in magnitudes.
        keywords['{}mag'.format(filter)] = row["{}mag".format(filter)]
        keywords['e_{}mag'.format(filter)] = row["e_{}mag".format(filter)]
    if not getattr(row['pmRA'],'mask',False):
        # Collect proper motions, if applicable.
        keywords['pmra'] = row['pmRA'] * u.mas / u.year
    if not getattr(row['pmDE'],'mask',False):
        # Collect proper motions, if applicable.
        keywords['pmdec'] = row['pmDE'] * u.mas / u.year
    for key in ['SpType', 'Name']:
        # Collect additional ASCII-formatted keys.
        keywords[key] = row[key].decode('ascii')
        
    # Build and return a target object.
    return Target(name=row['SimbadName'].decode('ascii'), position=position, _keywords=keywords)

def main():
    """Main function"""
    import argparse
    parser = argparse.ArgumentParser(description="Create a starlist of UKIRT Faint Stanadards")
    parser.add_argument("-o", "--output", type=argparse.FileType("w"), default="-", help="Output file")
    parser.add_argument("--fs", action='store_false', dest='full', help="Collect FS stars only")
    opt = parser.parse_args()
    date = datetime.datetime.now().isoformat()
    
    try:
        from astroquery.vizier import Vizier
    except ImportError as e:
        parser.error("Requires 'astroquery' to be installed.\n{0!r}".format(e))
    table = Vizier.get_catalogs([VIZIER_UKIRT_FS_CATALOG])[0]
    
    tl = TargetList([row_to_target(row) for row in table if (str(row['SimbadName'].decode('ascii')).startswith("FS") or opt.full)])
    tl.sort(key = lambda t : t.position.ra)
    
    # Write a header for the starlist so that we know the source.
    opt.output.write("# UKIRT Faint Standard Stars\n")
    opt.output.write("# Data from VIZIER catalog '{0:s}'\n".format(VIZIER_UKIRT_FS_CATALOG))
    opt.output.write("# ADS Reference: '{0:s}'\n".format(VIZIER_ADS_REFERENCE))
    opt.output.write("# Data collected on {0:s}\n".format(date))
    tl.to_starlist(opt.output)

if __name__ == '__main__':
    main()