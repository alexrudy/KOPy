# -*- coding: utf-8 -*-

from .base import Region

import astropy.units as u

__all__ = ['Box', 'Circle', 'Ruler']

class Box(Region):
    """A box shape."""
    
    TEMPLATE = 'box({ra:s},{dec:s},{width.value:f}",{height.value:f}",{position_angle.value:f})'
    
    def __init__(self, position, width, height, position_angle, **kwargs):
        super(Box, self).__init__(position, extra_kw=kwargs)
        self.width = u.Quantity(width, u.arcsec)
        self.height = u.Quantity(height, u.arcsec)
        self.position_angle = u.Quantity(position_angle, u.degree)
        
    
class Circle(Region):
    """A circle shape"""
    def __init__(self, position, radius, **kwargs):
        super(Circle, self).__init__(position, extra_kw=kwargs)
        self.radius = radius
        
    TEMPLATE = 'circle({ra:s},{dec:s},{radius.value:.2f}")'
    
class Ruler(Region):
    """A region which makes a ruler."""
    def __init__(self, start, stop, **kwargs):
        super(Ruler, self).__init__(start, extra_kw=kwargs)
        self.start = start
        self.stop = stop
        
    TEMPLATE = "# ruler({start_ra:s},{start_dec:s},{end_ra:s},{end_dec:s}) ruler=fk5 arcsec"
    
    def get_parameters(self):
        """Get the parameters for this object."""
        data = super(Ruler, self).get_parameters()
        data['start_ra'] = self.ra
        data['start_dec'] = self.dec
        data['end_ra'] = self.stop.ra.to_string(u.hourangle, sep=":", pad=False)
        data['end_dec'] = self.stop.dec.to_string(u.degree, sep=":", pad=False)
        return data
    

