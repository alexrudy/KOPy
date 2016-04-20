# -*- coding: utf-8 -*-

import numpy as np

import astropy.units as u

from astropy.coordinates import ICRS, SkyCoord
from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates.baseframe import (BaseCoordinateFrame, UnitSphericalRepresentation, CartesianRepresentation,
    FrameAttribute, QuantityFrameAttribute, RepresentationMapping, frame_transform_graph)
from astropy.coordinates.angles import rotation_matrix

class CoordinateLocationAttribute(FrameAttribute):
    """A frame attribute which is a coordinates object."""
    
    def convert_input(self, value):
        """
        Checks that the input is a SkyCoord with the necessary units (or the
        special value ``None``).

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        if value is None:
            return None, False
        elif isinstance(value, ICRS):
            return value, False
        else:
            if not hasattr(value, 'transform_to'):
                raise ValueError('"{0}" was passed into an '
                                 'CoordinateLocationAttribute, but it does not have '
                                 '"transform_to" method'.format(value))
            icrsobj = value.transform_to(ICRS)
            if isinstance(icrsobj, SkyCoord):
                return icrsobj.frame, True
            return icrsobj, True
    
    

class AstrometricFrame(BaseCoordinateFrame):
    """
    An instrument-relative frame, which has instrument X and Y

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    X : `Angle`, optional, must be keyword
    Y : `Angle`, optional, must be keyword
    positionangle : `Quantity`, optional, must be keyword
        The position angle of the coordinate frame.
    origin : `str`, optional, must be one of "imag" or "spec".
        Specificies the imager or spectograph pointing origin.

    """
    default_representation = UnitSphericalRepresentation

    frame_specific_representation_info = {
        'unitspherical': [RepresentationMapping('lat', 'Y'),
                          RepresentationMapping('lon', 'X')]
    }
    
    origin = CoordinateLocationAttribute(default=None)
    
    def at_origin(self, origin):
        """Return this frame at a new origin"""
        attrs = {}
        for name, value in self.get_frame_attr_names().items():
            attrs[name] = getattr(self, name, value)
        attrs['origin'] = origin
        return self.__class__(**attrs)

@frame_transform_graph.transform(FunctionTransform, AstrometricFrame, AstrometricFrame)
def astrometric_to_astrometric(from_astrometric_coord, to_astrometric_frame):
    """Transform between two osiris frames."""
    
    # If both frames have on-sky positions, then the transform should happen relative to both origins.
    if (from_astrometric_coord.origin is not None) and (to_astrometric_frame.origin is not None):
        return from_astrometric_coord.transform_to(ICRS).transform_to(to_astrometric_frame)
    
    # Otherwise, the transform occurs just by setting the new origin.
    return to_astrometric_frame.realize_frame(from_astrometric_coord.cartesian)
    
    
@frame_transform_graph.transform(FunctionTransform, ICRS, AstrometricFrame)
def icrs_to_astrometric(icrs_coord, astrometric_frame):
    """Convert an ICRS coordinate to an Astrometric frame."""
        
    # Define rotation matricies along the position angle vector, and
    # relative to the origin.
    mat2 = rotation_matrix(-astrometric_frame.origin.dec, 'y')
    mat3 = rotation_matrix(astrometric_frame.origin.ra, 'z')
    R = mat2 * mat3
    
    xyz = icrs_coord.cartesian.xyz
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)
    representation = CartesianRepresentation(xyz)
    return astrometric_frame.realize_frame(representation)
    
@frame_transform_graph.transform(FunctionTransform, AstrometricFrame, ICRS)
def astrometric_to_icrs(astrometric_coord, icrs_frame):
    """Convert an OSIRIS frame coordinaate to an ICRS"""
    
    # Define rotation matricies along the position angle vector, and
    # relative to the origin.
    mat2 = rotation_matrix(-astrometric_coord.origin.dec, 'y')
    mat3 = rotation_matrix(astrometric_coord.origin.ra, 'z')
    R = mat2 * mat3
    
    Rinv = np.linalg.inv(R)
    xyz = astrometric_coord.cartesian.xyz
    orig_shape = xyz.shape
    xyz = Rinv.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)
    
    representation = CartesianRepresentation(xyz)
    return icrs_frame.realize_frame(representation)