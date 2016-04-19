# -*- coding: utf-8 -*-
"""OSIRIS custom regions."""

import astropy.units as u
from astropy.coordinates import SkyCoord, BaseCoordinateFrame, ICRS

from ...regions import Box, Circle, Region, Ruler, RegionFile
from .coords import OsirisInstrumentFrame

def get_instrument_position_angle(coord):
    """Get the instrument position angle."""
    if isinstance(coord, SkyCoord):
        frame = coord.frame
    elif isinstance(coord, BaseCoordinateFrame):
        frame = coord
    else:
        raise TypeError("Must pass a coordinate frame.")
        
    if 'instrument_position_angle' in frame.get_frame_attr_names():
        return coord.instrument_position_angle
    else:
        return 0 * u.degree

class ImagerFOV(Box):
    """Imager FOV"""
    
    IMAGER_PA = 45 * u.degree
    IMAGER_SIZE = 20.4 * u.arcsec
    
    def __init__(self, position, **kwargs):
        position_angle = self.IMAGER_PA + get_instrument_position_angle(position)
        super(ImagerFOV, self).__init__(position, position_angle=position_angle,
            width=self.IMAGER_SIZE, height=self.IMAGER_SIZE, **kwargs)
    

def get_spec_FOV(filter, scale):
    """Compute the spectrograph field of view."""
    if filter.endswith("bb"):
        (w, h) = (16, 64)
    else:
        raise ValueError("I haven't comptued these scales yet...")
    return (w * scale, h * scale)

class SpecFOV(Box):
    """Spectrograph FOV"""
    
    def __init__(self, position, filter, scale, **kwargs):
        width, height = get_spec_FOV(filter, scale)
        position_angle = get_instrument_position_angle(position)
        super(SpecFOV, self).__init__(position, position_angle=position_angle,
            width=0.56 * u.arcsec, height=2.24*u.arcsec, **kwargs)
            
        
    

def create_region_from_DDF(ddf, target, guidestar=None):
    """Create and return a region file from a DDF."""
    target_frame = ddf.dataset.dithers.frame.at_origin(target)
    imager_frame = target_frame.at_origin(target_frame.imager)
    spec_frame   = target_frame.at_origin(target_frame.spectrograph)
    
    regions = RegionFile()
    regions.append(Circle(target, 0.5*u.arcsec, text=ddf.dataset.object.name, color="red"))
    if guidestar is not None:
        regions.append(Circle(guidestar, 0.5*u.arcsec, text="GS", color='white'))
        regions.append(Ruler(target, guidestar, color="white"))
    
    for i, dither in enumerate(ddf.dataset.dithers.positions):
        imager_pos = dither.position.transform_to(imager_frame)
        spec_pos = dither.position.transform_to(spec_frame)
        color = "cyan" if dither.sky else "green"
        if ddf.dataset.imag.enabled:
            regions.append(ImagerFOV(imager_pos, text="Imager {:d}".format(i+1), color=color))
        regions.append(SpecFOV(spec_pos, filter=ddf.dataset.spec.filter, 
            scale=ddf.dataset.spec.scale, text="Spec {:d}".format(i+1), color=color))
    
    size = 60 * u.arcsec if ddf.dataset.laser else 30 * u.arcsec
    position_angle = target_frame.instrument_position_angle
    if target_frame.pointing_origin == "spec":
        position_angle += 45 * u.deg
    regions.append(Box(target, size*2, size*2, position_angle = position_angle, color="white"))
    
    return regions
    