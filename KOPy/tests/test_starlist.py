import pytest
from .. import starlist
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np
from astropy.tests.helper import assert_quantity_allclose

def test_starlist_parse_line(starlist_line):
    """Ensure that we can parse a starlist line."""
    name, pos, kw = starlist.parse_starlist_line(starlist_line)
    assert kw['lgs'] == '1'
    assert kw['skip'] == '3'
    assert name.replace("_", " ") == "PG 0026+129"

def test_starlist_verify_line(starlist_line):
    """Check verification"""
    m = starlist.verify_starlist_line(starlist_line)
    assert not len(m)
    
def test_starlist_leading_whitespace():
    """Test starlist with leading whitespace."""
    starlist_line = "        tt020    09 00 20.4470  39 04 03.660 2000 b-r=0.9 b-v=0.5 rmag=13.66 sep=58.40 pa=265"
    name, pos, keywords = starlist.parse_starlist_line(starlist_line)
    assert name.strip() == "tt020"
    assert name == "        tt020"
    np.testing.assert_allclose(pos.ra.radian, Angle("09 00 20.4470", unit=u.hourangle).radian)
    np.testing.assert_allclose(pos.dec.radian, Angle("39 04 03.660", unit=u.degree).radian)
    
def test_starlist_read_from_fn(starlist_filename):
    """Test reading a starlist from a filename"""
    starlist.parse_starlist(starlist_filename)
    
def test_starlist_read_from_fd(starlist_fd):
    """Test reading a starlist from an open file"""
    starlist.parse_starlist(starlist_fd)
    
    
@pytest.mark.parametrize("keyword,value,units",[
    ("Jmag=10.2", 10.2 * u.mag, u.mag),
    ("B-Vmag=0.04", 0.04 * u.mag, u.mag),
    ("pmra=0.01", 0.01 * (u.hourangle / 3600 / u.year), (u.hourangle / 3600 / u.year)),
    ("pmdec=0.02", 0.02 * (u.arcsec / u.year), (u.arcsec / u.year)),
    ("dra=1.1", 1.1 * (u.hourangle / 3600 / u.hr), (u.hourangle / 3600 / u.hr)),
    ("ddec=2.0", 2.0 * (u.arcsec / u.hr), (u.arcsec / u.hr)),
    ("rotdest=70", 70 * u.degree, u.degree),
    ("raoffset=10", 10 * u.arcsec, u.arcsec),
    ("decoffset=10", 10 * u.arcsec, u.arcsec),
])
def test_starlist_parse_keywords(keyword, value, units, starlist_line):
    """Test that quantities are parsed correctly in keywords."""
    starlist_line += " " + keyword
    kwname = keyword.split("=",1)[0]
    name, position, keywords = starlist.parse_starlist_line(starlist_line)
    assert kwname in keywords
    assert keywords[kwname].unit.is_equivalent(units)
    assert_quantity_allclose(keywords[kwname], value)