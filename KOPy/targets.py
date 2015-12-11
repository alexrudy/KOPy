# -*- coding: utf-8 -*-
#
#  targets.py
#  KOPy
#
#  Created by Alexander Rudy on 2015-08-07.
#  Copyright 2015 Alexander Rudy. All rights reserved.
#
"""
:mod:`targets` provides an object-oriented way to work with starlists in the Keck starlist format.

For example::
    
    >>> t = Target.from_starlist("S2CM006571      05 04 37.23 -19 37 04.58 2000 b-r=2.38 b-v=0.57 rmag=14.73")
    >>> t
    <Target 'S2CM006571'@'05h04m37.23s -19d37m04.58s' rmag=14.73, b-r=2.38, b-v=0.57>
    >>> t.to_starlist()
    'S2CM006571      05 04 37.230 -19 37 04.580 2000 rmag=14.73 b-r=2.38 b-v=0.57'
    >>> tl = TargetList([t])
    >>> tl
    [<Target 'S2CM006571'@'05h04m37.23s -19d37m04.58s' rmag=14.73, b-r=2.38, b-v=0.57>]
    >>> tl.to_starlist(filename=None)
    'S2CM006571      05 04 37.230 -19 37 04.580 2000 rmag=14.73 b-r=2.38 b-v=0.57\\n'

"""

import six
import astropy.units as u
import collections
import io
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, Column, MaskedColumn
from .starlist import (parse_starlist, parse_starlist_line, 
    format_starlist_line, format_keywords, format_starlist_position)

__all__ = ['Target', 'TargetList']

class Target(object):
    """A single target object, with a position, name, and keyword arguments.
    
    Parameters
    ----------
    name : string
        Name of the target.
    position : :class:`~astropy.coordinates.SkyCoord`
        Sky coordinate position.
    keywords : Key-value pairs
        Additional keywords can be passed to the constructor and will be used
        to populate the keyword dictionary. To preserve the order of keywords
        passed to the constructor, pass an ordered dict or similar to the
        ``_keywords`` argument.
    
    Examples
    --------
    
    Construct a target from a name resolved via SIMBAD:
        
        >>> Target.from_name("M1") # doctest: +SKIP
        <Target 'M1'@...>
        >>> t = Target("Galaxy", SkyCoord(10, 12, unit=(u.deg, u.deg), frame='icrs'))
        >>> t
        <Target 'Galaxy'@...>
        >>> t.name
        'Galaxy'
        >>> t2 = Target("GuideStar", t.position, rmag=12)
        >>> t2.keywords['vmag'] = 15
        >>> t2.keywords['rmag']
        12
        >>> t2
        <Target 'GuideStar'@... rmag=12.00, vmag=15.00>
    
    Using a constructed target, make a starlist item:
        
        >>> t2.to_starlist()
        'GuideStar       00 40 00.002 +12 00 00.006 2000 rmag=12.00 vmag=15.00'
    
    """
    
    #TODO: Can't really use slots here. Does it matter?
    # It breaks pickling and more :(
    # __slots__ = ('name', 'position', 'keywords')
    
    name = None
    """Name of the target"""

    position = None
    """Position of the target, an :class:`~astropy.coordinates.SkyCoord` object."""

    keywords = {}
    """Keyword-value pairs from the starlist, as a dictionary."""
    
    def __init__(self, name, position, _keywords=dict(), **kwargs):
        super(Target, self).__init__()
        self.name = str(name)
        self.position = position if isinstance(position, SkyCoord) else SkyCoord(position)
        self.keywords = collections.OrderedDict()
        self.keywords.update(_keywords)
        self.keywords.update(kwargs)
        
        
    def __getattr__(self, key):
        """Delegate attributes to keywords when necessary."""
        if hasattr(self, 'keywords') and key in self.keywords:
            return self.keywords[key]
        raise AttributeError("'{:s}' has no attribute '{:s}'".format(self.__class__.__name__, key))
        
    def __setattr__(self, key, value):
        """Set an attribute as a keyword."""
        if not hasattr(self, key):
            self.keywords[key] = value
        else:
            super(Target, self).__setattr__(key, value)
        
    def __repr__(self):
        """Represent a target."""
        return "<{0:s} '{1:s}'@'{2:s}' {3:s}>".format(self.__class__.__name__, 
            self.name, self.pos_string(), ", ".join(self._repr_keywords_()))
        
    def __eq__(self, other):
        """Equality."""
        return (self.name == other.name) and (self.position == other.position)
        
    def _repr_keywords_(self):
        """Clean representation of keywords."""
        return format_keywords(self.keywords)
        
    def pos_string(self):
        """Position string, in a nicely formatted way."""
        return self.position.to_string('hmsdms')
        
    def to_starlist(self, **kwargs):
        """Return a starlist line."""
        return format_starlist_line(self.name, self.position, self.keywords, **kwargs)
        
    @classmethod
    def from_starlist(cls, line):
        """Parse a single line from a starlist into the Target data structure."""
        name, position, kw = parse_starlist_line(line)
        return cls(name=name, position=position, _keywords=kw)
        
    @classmethod
    def from_name(cls, name):
        """Use the SkyCoord name resolution (which works via Simbad) to find the coordinates of a name."""
        position = SkyCoord.from_name(name)
        return cls(name=name, position=position)
    

class TargetList(collections.MutableSequence):
    """A target list.
    
    Target lists are lists of :class:`Target` objects, with a few special methods.
    
    Parameters
    ----------
    iterable : 
    
    """
    def __init__(self, iterable=None):
        super(TargetList, self).__init__()
        self.__data = []
        if iterable is not None:
            self.extend(iterable)
    
    def __repr__(self):
        """Represent the target list."""
        return repr(self.__data)
    
    @classmethod
    def _type_check(cls, value):
        """Type check the value."""
        if not isinstance(value, Target):
            raise TypeError("{0:s} must contain only subclasses of {1:s}".format(
                cls.__name__, Target.__name__
            ))
        return value
    
    def __setitem__(self, key, value):
        """Ensure type consistency!"""
        return self.__data.__setitem__(key, self._type_check(value))
    
    if six.PY2:
        def __getslice__(self, start, end):
            """Override the slice operator in python2."""
            return self.__getitem__(slice(start, end))
    
    def __getitem__(self, key):
        """Ensure that we always get back a TargetList"""
        
        # Support indexing by name.
        if isinstance(key, six.string_types):
            for t in self:
                if t.name == key:
                    return t
            else:
                raise KeyError("No target with name '{0:s}' found".format(key))
        
        # Support regular list indexing, but turn it into a TargetList
        # if we would otherwise return a list.
        r = self.__data.__getitem__(key)
        if isinstance(r, list):
            return self.__class__(r)
        return r
        
    def __delitem__(self, key):
        """Delete an item by key."""
        if isinstance(key, six.string_types):
            for i, t in self:
                if t.name == key:
                    key = i
                    break
            else:
                raise ValueError("No target with name '{0:s}' found".format(key))
        return self.__data.__delitem__(key)
        
    def __add__(self, item):
        """Add two TargetList objects together."""
        return self.__class__(self.__data.__add__(item))
        
    def __mul__(self, value):
        """Multiply"""
        return self.__class__(self.__data.__mul__(value))
    
    def __rmul__(self, value):
        """Reverse multiply."""
        return self.__class__(self.__data.__rmul__(value))
        
    def __len__(self):
        """Length, from the underlying list."""
        return self.__data.__len__()
        
    def sort(self, *args, **kwargs):
        """Sort the list."""
        return self.__data.sort(*args, **kwargs)
        
    def insert(self, index, item):
        """Insert an item, and check type."""
        return self.__data.insert(index, self._type_check(item))
    
    @classmethod
    def from_starlist(cls, filename):
        """From a starlist"""
        return cls(Target(name, position, _keywords=kw) for name, position, kw in parse_starlist(filename))
    
    @classmethod
    def from_table(cls, table):
        """Reconstruct the starlist from a table created by self.table()"""
        if "Name" in table.colnames:
            names = table['Name']
            table.remove_column('Name')
        else:
            raise ValueError("Table must have a 'Name' column.")
        if "Position" in table.colnames:
            positions = SkyCoord(table['Position'])
            table.remove_column('Position')
        elif "RA" in table.colnames and "Dec" in table.colnames:
            positions = SkyCoord(table['RA'], table['Dec'])
            table.remove_columns(['RA', 'Dec'])
        else:
            raise ValueError("Table must have either 'Position' or 'RA' and 'Dec' columns.")
        return cls(Target(name, position, _keywords=[(k,row[k]) for k in row.colnames if row[k] is not np.ma.masked ]) for name, position, row in zip(names, positions, table))
    
    @property
    def names(self):
        """Target names."""
        return [ t.name for t in self ]
        
    def catalog(self):
        """Make a single SkyCoord object for all targets."""
        return SkyCoord([t.position.transform_to('icrs') for t in self])
        
    def table(self, coord_mixin=False):
        """Create a table object which represents this target list."""
        
        cols = collections.OrderedDict()
        def _make_row(t):
            d = collections.OrderedDict(Name=t.name)
            d.update(t.keywords)
            cols.update((k, None) for k in d.keys())
            return d
        
        
        data = [ _make_row(t) for t in self ]
        columns = []
        for c in cols:
            ma = np.ma.MaskedArray([d.get(c, np.ma.masked) for d in data])
            columns.append(MaskedColumn(ma, name=c, mask=ma.mask))
        
        t = Table(columns, masked=True)
        if coord_mixin:
            t['Position'] = self.catalog()
            t['Position'].format = lambda c : c.to_string('hmsdms')
        else:
            catalog = self.catalog()
            t.add_column(MaskedColumn(catalog.ra.to(u.hourangle), name='RA', format=lambda c : Angle(c).to_string(), mask=np.zeros(catalog.shape, dtype=np.bool)), index=1)
            t.add_column(MaskedColumn(catalog.dec, name="Dec", format=lambda c : Angle(c).to_string(), mask=np.zeros(catalog.shape, dtype=np.bool)), index=2)
        return t
        
    def _to_starlist_stream(self, file, **kwargs):
        """Write to a stream"""
        for t in self:
            file.write(t.to_starlist(**kwargs))
            file.write("\n")
        
    def to_starlist(self, filename, mode='w', **kwargs):
        """Write to a starlist file."""
        if filename is None:
            return "\n".join([t.to_starlist(**kwargs) for t in self]) + "\n"
        if isinstance(filename, io.IOBase) or hasattr(filename, 'write'):
            stream = filename
            self._to_starlist_stream(stream, **kwargs)
        else:
            with open(filename, mode) as stream:
                self._to_starlist_stream(stream, **kwargs)
        