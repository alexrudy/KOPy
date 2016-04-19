import pytest

import os
import pickle

from ..targets import Target, TargetList
from astropy.tests.helper import assert_quantity_allclose
from astropy.tests.helper import pickle_protocol
import astropy.units as u
from astropy.coordinates import SkyCoord
from six.moves import StringIO

def pickle_roundtrip(original, protocol):
    """Round trip something via Pickle"""
    f = pickle.dumps(original, protocol=protocol)
    unpickled = pickle.loads(f)
    return unpickled
    
def assert_coord_allclose(actual, desired, rtol=1.e-7, atol=None):
    """Check that coordinates are allclose."""
    actual = SkyCoord(actual)
    desired = SkyCoord(desired)
    assert_quantity_allclose(actual.ra, desired.ra, rtol=rtol, atol=atol)
    assert_quantity_allclose(actual.dec, desired.dec, rtol=rtol, atol=atol)
    assert_quantity_allclose(actual.distance, desired.distance, rtol=rtol, atol=atol)
    
def assert_target_allclose(actual, desired, rtol=1.e-7, atol=None):
    """docstring for assert_target_allclose"""
    assert actual.name == desired.name
    for keyword in actual.keywords:
        assert actual.keywords[keyword] == desired.keywords[keyword]
        assert getattr(actual, keyword) == getattr(desired, keyword)
    assert_coord_allclose(actual.position, desired.position)

@pytest.fixture(scope='module')
def targetlist():
    """An example target list, quite simple"""
    import pkg_resources
    return TargetList.from_starlist(pkg_resources.resource_filename(__name__, 'data/small_starlist.txt'))

@pytest.fixture
def target():
    """A target option."""
    return Target(name="MyTarget", position="1h4m3s +5d4m9s", param="value")

def test_read_starlist(starlist_filename):
    """Read a starlist filename."""
    tl = TargetList.from_starlist(starlist_filename)
    assert len(tl)
    
def test_index_starlist(targetlist):
    """Test indexing into a starlist."""
    assert len(targetlist) == 15
    
    s = targetlist[1:3]
    assert len(s) == 2
    assert isinstance(s, TargetList)
    
    s = targetlist[3:1:-1]
    assert len(s) == 2
    assert isinstance(s, TargetList)
    
    t = targetlist[1]
    assert isinstance(t, Target)
    
    t = targetlist['ring neb']
    assert t.name == 'ring neb'
    
    with pytest.raises(KeyError):
        targetlist['OTHERNAME']
    
def test_setitem_starlist(targetlist, target):
    """Test target item starlist."""
    targetlist[1] = target
    assert targetlist[1].name == target.name
    
def test_starlist_type_enforcement():
    """Test type enforcement."""
    
    with pytest.raises(TypeError):
        tl = TargetList([1])
        
    tl = TargetList()
    with pytest.raises(TypeError):
        tl.append(1)
        
    tl = TargetList()
    with pytest.raises(TypeError):
        tl.extend([1,2,3])
        
def test_starlist_add(targetlist, target):
    """Test addition of starlists."""
    tl = targetlist + []
    assert isinstance(tl, TargetList)
    
    tl = targetlist + [target]
    assert isinstance(tl, TargetList)
    
    with pytest.raises(TypeError):
        tl = targetlist + [1]
        
def test_starlist_mul(target):
    """Test starlist multiplication."""
    tl = TargetList([target]) * 2
    assert isinstance(tl, TargetList)
    assert len(tl) == 2
    
    tl = 4 * TargetList([target])
    assert isinstance(tl, TargetList)
    assert len(tl) == 4
    
def test_starlist_names(target):
    """Test target names."""
    tl = TargetList([target, target])
    assert tl.names == [target.name, target.name]
    
def test_starlist_catalog(target):
    """Test starlist catalog"""
    tl = TargetList([target, target])
    assert_coord_allclose(tl.catalog(), SkyCoord([target.position, target.position]))
    
def test_starlist_leading_whitespace(tricky_tt):
    """Test starlist with leading whitespace."""
    starlist_line, pos = tricky_tt
    t = Target.from_starlist(starlist_line)
    assert t.name.strip() == "tt020"
    assert t.name == "        tt020"
    assert_coord_allclose(t.position, pos)
    assert t.to_starlist().startswith("        tt020")
    
@pytest.mark.xfail
def test_targetlist_table_roundtrip(targetlist):
    """Round trip via table."""
    tl = TargetList.from_table(targetlist.table())
    assert tl.names == targetlist.names
    
def test_targetlist_starlist_roundtrip(targetlist, tmpdir):
    """Target list round trip via starlist."""
    targetlist.to_starlist(str(tmpdir.join('starlist.txt')))
    tl = TargetList.from_starlist(str(tmpdir.join('starlist.txt')))
    assert tl.names == targetlist.names
    
def test_targetlist_pickle_roundtrip(target, pickle_protocol):
    """Test for pickleing round-trip."""
    unpickled = pickle_roundtrip(target, pickle_protocol)
    assert_target_allclose(unpickled, target)