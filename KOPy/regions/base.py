# -*- coding: utf-8 -*-

import astropy.units as u

HEADER_TEMPATE = """
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
"""[1:-1]

class RegionFile(list):
    """Region File"""
    def __init__(self, *items):
        super(RegionFile, self).__init__(items)
        
    def __call__(self):
        """Render the region file."""
        for line in HEADER_TEMPATE.splitlines():
            yield line
        for item in self:
            yield item()

class Region(object):
    """A region."""
    
    def __init__(self, position, extra_kw):
        super(Region, self).__init__()
        self.position = position.transform_to("fk5")
        self._extra = extra_kw
    
    @property
    def ra(self):
        """RA"""
        return self.position.ra.to_string(u.hourangle, sep=":", pad=False)
    
    @property
    def dec(self):
        """Declination"""
        return self.position.dec.to_string(u.degree, sep=":", pad=False)
    
    TEMPLATE = ""
    
    def get_parameters(self):
        """Get the parameters for the template string."""
        data = dict()
        data.update(vars(self))
        data['ra'] = self.ra
        data['dec'] = self.dec
        return data
    
    def get_extras(self):
        """Get the extras."""
        extra = []
        for keyword, value in self._extra.items():
            if keyword == "text":
                extra.append("{0}={{{1}}}".format(keyword, value))
            else:
                extra.append("{0}={1}".format(keyword, value))
        return " ".join(extra)
    
    def __call__(self):
        """Render the region."""
        core = self.TEMPLATE.format(**self.get_parameters())
        extra = self.get_extras()
        if "#" not in core:
            core += " # "
        return core + " " + extra
        
    

        