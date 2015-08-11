import pytest
from .. import starlist

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
    
def test_starlist_read_from_fn(starlist_filename):
    """Test reading a starlist from a filename"""
    starlist.parse_starlist(starlist_filename)
    
def test_starlist_read_from_fd(starlist_fd):
    """Test reading a starlist from an open file"""
    starlist.parse_starlist(starlist_fd)
    