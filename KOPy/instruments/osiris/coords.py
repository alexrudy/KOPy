# -*- coding: utf-8 -*-

import numpy as np

import astropy.units as u

from astropy.coordinates import ICRS, SkyCoord
from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates.baseframe import (BaseCoordinateFrame, UnitSphericalRepresentation, CartesianRepresentation,
    FrameAttribute, QuantityFrameAttribute, RepresentationMapping, frame_transform_graph, SphericalRepresentation)
from astropy.coordinates.angles import rotation_matrix

from ..coords import AstrometricFrame, CoordinateLocationAttribute

__all__ = ['OsirisInstrumentFrame']

_OSIRIS_IMAGER_OFFSET = 19.4 * u.arcsec
_OSIRIS_IMAGER_RELATIVE_PA = 135.0 * u.degree

class OsirisInstrumentFrame(AstrometricFrame):
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
    
    pointing_origin = FrameAttribute(default="spec")
    rotation = QuantityFrameAttribute(default=0 * u.radian, unit=u.radian)
    
    @property
    def imager(self):
        """The origin of the imager in this frame."""
        if self.pointing_origin == "imag":
            return self.realize_frame(UnitSphericalRepresentation(0.0 * u.deg, 0.0 * u.deg))
        else:
            xyz = UnitSphericalRepresentation(_OSIRIS_IMAGER_OFFSET, 0.0 * u.deg).to_cartesian().xyz
            m1 = rotation_matrix(_OSIRIS_IMAGER_RELATIVE_PA, 'x')
            R = m1
            xyz = R.dot(xyz)
            return self.realize_frame(CartesianRepresentation(*xyz))
        
    @property
    def spectrograph(self):
        """The origin of the spectrograph in this frame."""
        if self.pointing_origin == "spec":
            return self.realize_frame(UnitSphericalRepresentation(0.0 * u.deg, 0.0 * u.deg))
        else:
            xyz = UnitSphericalRepresentation(0.0 * u.deg, -_OSIRIS_IMAGER_OFFSET).to_cartesian().xyz
            R = rotation_matrix(_OSIRIS_IMAGER_RELATIVE_PA + self.instrument_position_angle - self.rotation, 'x')
            xyz = R.dot(xyz)
            return self.realize_frame(CartesianRepresentation(*xyz))
        

@frame_transform_graph.transform(FunctionTransform, OsirisInstrumentFrame, OsirisInstrumentFrame)
def osiris_to_osiris(from_osiris_coord, to_osiris_frame):
    """Transform between two osiris frames."""
    
    # If both frames have on-sky positions, then the transform should happen relative to both origins.
    if (from_osiris_coord.origin is not None) and (to_osiris_frame.origin is not None):
        return from_osiris_coord.transform_to(ICRS).transform_to(to_osiris_frame)
    
    # Otherwise, the transform occurs just by setting the new origin and rotating
    xyz = from_osiris_coord.cartesian.xyz
    R = rotation_matrix(-to_osiris_frame.rotation + from_osiris_coord.rotation, 'x')
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)
    
    representation = CartesianRepresentation(xyz)
    return to_osiris_frame.realize_frame(representation)
    
    
@frame_transform_graph.transform(FunctionTransform, ICRS, OsirisInstrumentFrame)
def icrs_to_osiris(icrs_coord, osiris_frame):
    """Convert an ICRS coordinate to an Astrometric frame."""
        
    # Define rotation matricies along the position angle vector, and
    # relative to the origin.
    mat0 = np.identity(3)
    mat0[1,1] = -1.0
    mat0 = np.matrix(mat0)
    mat1 = rotation_matrix(-osiris_frame.rotation, 'x')
    mat2 = rotation_matrix(-osiris_frame.origin.dec, 'y')
    mat3 = rotation_matrix(osiris_frame.origin.ra, 'z')
    R = mat0 * mat1 * mat2 * mat3
    
    xyz = icrs_coord.cartesian.xyz
    orig_shape = xyz.shape
    xyz = R.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)
    # xyz[1] *= -1
    representation = CartesianRepresentation(xyz)
    return osiris_frame.realize_frame(representation)
    
@frame_transform_graph.transform(FunctionTransform, OsirisInstrumentFrame, ICRS)
def osiris_to_icrs(osiris_coord, icrs_frame):
    """Convert an OSIRIS frame coordinaate to an ICRS"""
    
    # Define rotation matricies along the position angle vector, and
    # relative to the origin.
    mat0 = np.identity(3)
    mat0[1,1] = -1.0
    mat0 = np.matrix(mat0)
    mat1 = rotation_matrix(-osiris_coord.rotation, 'x')
    mat2 = rotation_matrix(-osiris_coord.origin.dec, 'y')
    mat3 = rotation_matrix(osiris_coord.origin.ra, 'z')
    R = mat0 * mat1 * mat2 * mat3
    
    Rinv = np.linalg.inv(R)
    xyz = osiris_coord.cartesian.xyz
    # xyz[1] *= -1
    orig_shape = xyz.shape
    xyz = Rinv.dot(xyz.reshape(xyz.shape[0], np.prod(xyz.shape[1:]))).reshape(orig_shape)
    
    representation = CartesianRepresentation(xyz)
    return icrs_frame.realize_frame(representation)