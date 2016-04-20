# -*- coding: utf-8 -*-
"""OSIRIS custom regions."""
import os
import argparse
import six
import glob
import subprocess

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
        
    if 'rotation' in frame.get_frame_attr_names():
        return coord.rotation
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
    scale = scale.to(u.arcsec/u.pix).value * u.arcsec
    return (w * scale, h * scale)

class SpecFOV(Box):
    """Spectrograph FOV"""
    
    def __init__(self, position, filter, scale, **kwargs):
        width, height = get_spec_FOV(filter, scale)
        position_angle = get_instrument_position_angle(position)
        super(SpecFOV, self).__init__(position, position_angle=position_angle,
            width=width, height=height, **kwargs)
            
        
    

def create_region_from_DDF(ddf, target, guidestar=None):
    """Create and return a region file from a DDF."""
    target_frame = OsirisInstrumentFrame(origin=target, rotation=ddf.dataset.dithers.position_angle)
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
    position_angle = target_frame.rotation
    if target_frame.pointing_origin == "spec":
        position_angle += 45 * u.deg
    regions.append(Box(target, size*2, size*2, position_angle = position_angle, color="white"))
    
    return regions
    

def ddf_open_ds9(target, filename, imdir="."):
    """Assemble the arguments to open DS9"""
    ds9args = ["ds9", '-view', 'layout', 'vertical']
    impath = os.path.join(os.path.abspath(imdir),"{}*.fits".format(target.name))
    images = glob.glob(impath)
    if not len(images):
        raise IOError("No images found in '{}'".format(impath))
    for idx, image in enumerate(images):
        if idx:
            ds9args += ['-frame', 'new']
        ds9args += [image, '-cmap', 'hsv', '-scale', 'log', '-bg', 'black']
    ds9args += ["-regions", "load", 'all', filename ]
    ds9args += ["-pan", "to",
        "{:s}".format(target.position.ra.to_string(u.hourangle, sep=":", pad=False)), 
        "{:s}".format(target.position.dec.to_string(u.degree, sep=":", pad=False)), "wcs", "fk5"]
    ds9args += ['-frame', 'match', 'wcs']
    ds9args += ['-frame', 'lock', 'wcs']
    ds9args += ['-frame', six.text_type(idx+1)]
    ds9args += ['-geometry', '1300x1100']
    print(" ".join(ds9args))
    subprocess.Popen(ds9args)

def main(args=None):
    """Main function for DDF manipulation."""
    parser = argparse.ArgumentParser(description="Produce a simple region file for help making OSIRIS DDFs")
    parser.add_argument('starlist', help='Starlist filename', type=six.text_type)
    parser.add_argument('ddf', help='DDF filename', type=six.text_type)
    parser.add_argument('target', help='Target name', type=six.text_type, nargs="?")
    
    parser.add_argument('--output', help='Output region name', type=argparse.FileType('w'))
    parser.add_argument('--ds9', help='Open with DS9', action='store_true')
    parser.add_argument('--imdir', help='Image directory', type=six.text_type, default=os.path.relpath(os.getcwd()))
    
    opt = parser.parse_args(args)
    
    from .ddf import DataDefinitionFile
    from ...targets import TargetList
    
    print("Parsing DDF '{}'".format(opt.ddf))
    ddf = DataDefinitionFile.from_file(opt.ddf)
    if opt.target is None:
        opt.target = ddf.dataset.name
    
    print("Parsing starlist '{}'".format(opt.starlist))
    targets = TargetList.from_starlist(opt.starlist)
    
    target = targets[opt.target]
    guidestar = targets[targets.index(target) + 1]
    if guidestar.position.separation(target.position) > 80 * u.arcsec:
        guidestar = None
    
    print("Creating region file.")
    regions = create_region_from_DDF(ddf, target.position, guidestar.position)
    
    if opt.output is None: 
        opt.output = open("{}-auto.reg".format(target.name), 'w')
    print("Saving region file to '{}'".format(opt.output.name))
    opt.output.write("\n".join(regions()))
    opt.output.flush()
    opt.output.close()
    
    if opt.ds9:
        print("Opening images in DS9 from {}".format(opt.imdir))
        try:
            ddf_open_ds9(target, opt.output.name, opt.imdir)
        except IOError as e:
            parser.error(e)
    