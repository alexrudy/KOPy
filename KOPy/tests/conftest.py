import pytest

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