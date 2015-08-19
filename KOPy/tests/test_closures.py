# -*- coding: utf-8 -*-
"""
Tests for the LGS closure parser.
"""

import pytest
from KOPy.closures import LCHRegions
from astropy.time import Time

@pytest.fixture
def closures_filename():
    """A filename for closures."""
    import pkg_resources
    return pkg_resources.resource_filename(__name__, "data/opensUnix150807.txt")
    

def test_parse_closures(closures_filename):
    """Test parse closures."""
    regions = LCHRegions.parse(closures_filename, date='2015-08-07')
    assert "lazer_zenith" in regions
    assert "eng341" in regions
    
    assert len(regions.keys())
    assert len(regions.values())
    
    lazer_zenith = regions["lazer_zenith"]
    assert lazer_zenith.open(Time('2015-08-07 11:50:29'))
    assert not lazer_zenith.open(Time('2015-08-07 11:50:32'))
    assert not lazer_zenith.open(Time('2015-08-07 11:50:30'))
    
    