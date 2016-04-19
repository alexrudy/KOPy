# -*- coding: utf-8 -*-
"""
Data Definition File
"""

import astropy.units as u
from astropy.coordinates import SkyCoord

from .coords import OsirisInstrumentFrame

class SpectrographParameters(object):
    """Spectrograph parameters"""
    def __init__(self, filter, scale, itime, coadds):
        super(SpectrographParameters, self).__init__()
        self.filter = filter
        self.scale = u.Quantity(scale, unit=u.arcsec/u.pix)
        self.itime = u.Quantity(itime, unit=u.s)
        self.coadds = int(coadds)
    
    @classmethod
    def parse(cls, xml):
        """Parse an XML tree for the dataset parameters."""
        scale = float(xml.get('scale').split('"')[0])
        return cls(filter=xml.get('filter'), scale=scale, itime=float(xml.get('itime')), 
            coadds=int(xml.get('coadds')))

class ImagerFrames(object):
    """An individual imager frame"""
    def __init__(self, filter, itime, coadds, repeats):
        super(ImagerFrames, self).__init__()
        self.filter = filter
        self.itime = u.Quantity(itime, unit=u.s)
        self.coadds = int(coadds)
        self.repeats = int(repeats)
    
    @classmethod
    def parse(cls, xml):
        """Parse an XML tree for the dataset parameters."""
        return cls(filter=xml.get('filter'), itime=float(xml.get('itime')), 
            coadds=int(xml.get('coadds')), repeats=int(xml.get('repeats')))

class ImagerParameters(object):
    """Imager parameters"""
    def __init__(self, frames, mode):
        super(ImagerParameters, self).__init__()
        self.frames = frames
        self.mode = mode
        
    @property
    def enabled(self):
        """Is the imager enabled."""
        return self.mode != "Disabled (Spec only)"
        
    @classmethod
    def parse(cls, xml):
        """Parse an XML tree for the dataset parameters."""
        frames = [ ImagerFrames.parse(frame) for frame in xml.findall('imagFrame') ]
        return cls(frames, mode=xml.get('mode'))

class DitherPosition(object):
    """A dither position"""
    def __init__(self, position, sky=False):
        super(DitherPosition, self).__init__()
        self.position = position
        self.sky = sky
        

class DitherPattern(object):
    """A collection of dither positions."""
    def __init__(self, frame, positions):
        super(DitherPattern, self).__init__()
        self.frame = frame
        self.positions = positions
        
    @property
    def imager(self):
        """Get the imager positions."""
        return SkyCoord([pos.position for pos in self.positions], frame=self.frame) + self.frame.imager
        
    @property
    def spectrograph(self):
        """Get the spectrograph positions."""
        return SkyCoord([pos.position for pos in self.positions], frame=self.frame) + self.frame.spectrograph
        
    @classmethod
    def parse(cls, xml):
        """Parse the XML for a series of dither positions."""
        coords = xml.get('coords')
        position_angle = float(xml.get('skyPA')) * u.degree
        units = u.Unit(xml.get('units'))
        frame = OsirisInstrumentFrame(pointing_origin="spec", instrument_position_angle=position_angle)
        positions = []
        for position in xml.findall('ditherPosition'):
            X = float(position.get('xOff')) * units
            Y = float(position.get('yOff')) * units
            sky = position.get('sky') == 'true'
            positions.append(DitherPosition(SkyCoord(X=X, Y=Y, frame=frame), sky=sky))
        return cls(frame, positions)

class DatasetParameters(object):
    """Dataset parameters"""
    def __init__(self, name, number, aomode, object, spec, imag, dithers):
        super(DatasetParameters, self).__init__()
        self.name = name
        self.number = int(number)
        self.aomode = aomode
        self.object = object
        self.spec = spec
        self.imag = imag
        self.dithers = dithers
        
        
    @property
    def laser(self):
        """Is this a laser mode?"""
        return "LGS" in self.aomode
    
    @classmethod
    def parse(cls, xml):
        """Parse an XML tree for the dataset parameters."""
        name = xml.get('name')
        number = int(xml.get('setnum'))
        aomode = xml.get('aomode')
        obj = Object.parse(xml.find('object'))
        spec = SpectrographParameters.parse(xml.find('spec'))
        imag = ImagerParameters.parse(xml.find('imag'))
        dither = DitherPattern.parse(xml.find('ditherPattern'))
        return cls(name, number, aomode, obj, spec, imag, dither)

class Object(object):
    """Object parameters"""
    def __init__(self, name):
        super(Object, self).__init__()
        self.name = name
        
    @classmethod
    def parse(cls, xml):
        """Parse an XML tree for the object parameter."""
        return cls("".join(xml.itertext()))

class DataDefinitionFile(object):
    """A container for the data definition file format."""
    def __init__(self, dataset):
        super(DataDefinitionFile, self).__init__()
        self.dataset = dataset
    
    @classmethod
    def from_file(cls, file):
        """Read the DDF from a file."""
        import xml.etree.ElementTree as ET
        return cls.parse(ET.parse(file))
        
    @classmethod
    def parse(cls, xml):
        """Parse the XML."""
        ds = DatasetParameters.parse(xml.find('dataset'))
        return cls(ds)
        