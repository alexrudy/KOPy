import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

starlist_tokens = [
    ["PG_0026+129     ","8 47 42.5 ","+34 45 04.3 ","2000 ","lgs=1 skip=3"],
    ["PG 0026+129     ","00 29 13.7000 ","13 16 03.720 ","2000 ","lgs=1 skip=3"],
    ["PG 0026+129     ","0 2 1 ","3 1 3 ","2000 ","lgs=1 skip=3"],
    ]

@pytest.fixture(params=starlist_tokens)
def tokens(request):
    """Starlist tokens"""
    return request.param
    
@pytest.fixture
def starlist_line(tokens):
    """Starlist token lines."""
    return "".join(tokens)


@pytest.fixture(params=[
    ("        tt020    09 00 20.4470  39 04 03.660 2000 b-r=0.9 b-v=0.5 rmag=13.66 sep=58.40 pa=265", SkyCoord("09 00 20.4470", "39 04 03.660", unit=(u.hourangle, u.degree))),
    ("        tt020    09 15 51.3420  44 19 59.860 2000 b-r=0.9 b-v=0.5 rmag=15.98 sep=36.47 pa=345", SkyCoord("09 15 51.3420", "44 19 59.860", unit=(u.hourangle, u.degree)))
])
def tricky_tt(request):
    """Return tricky TT lines to parse."""
    return request.param

@pytest.fixture(params=[
    pytest.mark.slow('data/big_starlist.txt'),
    'data/small_starlist.txt',
])
def _starlist_fn(request):
    """Starlist filename"""
    return request.param

@pytest.fixture
def starlist_filename(_starlist_fn):
    """The starlist file."""
    import pkg_resources
    return pkg_resources.resource_filename(__name__, _starlist_fn)

@pytest.fixture
def starlist_fd(_starlist_fn):
    """Starlist file descriptor"""
    import pkg_resources
    return pkg_resources.resource_stream(__name__, _starlist_fn)