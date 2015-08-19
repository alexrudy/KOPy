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
"""

import six
import astropy.units as u
import collections
import io
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table, Column, MaskedColumn
from .starlist import parse_starlist, parse_starlist_line

class Target(object):
    """A single target object, with a position, name, and keyword arguments."""
    
    keywords = {}
    
    def __init__(self, name, position, _keywords=dict(), **kwargs):
        super(Target, self).__init__()
        self.name = str(name)
        self.position = SkyCoord(position)
        self.keywords = collections.OrderedDict()
        self.keywords.update(_keywords)
        self.keywords.update(kwargs)
        
        
    def __getattr__(self, key):
        """Delegate attributes to keywords when necessary."""
        if key in self.keywords:
            return self.keywords[key]
        raise AttributeError("'{:s}' has no attribute '{:s}'".format(self.__class__.__name__, key))
        
    def __repr__(self):
        """Represent a target."""
        return "<{0:s} '{1:s}'@'{2:s}' {3:s}>".format(self.__class__.__name__, self.name, self.pos_string(), ", ".join(self._repr_keywords_()))
        
    def __eq__(self, other):
        """Equality."""
        return (self.name == other.name) and (self.position == other.position)
        
    def _repr_keywords_(self):
        """Clean representation of keywords."""
        if not len(self.keywords):
            return [""]
        keywords = []
        for key, value in self.keywords.items():
            if "mag" == key[-3:]:
                keywords.insert(0, "{key}={value:.2f}".format(key=key, value=float(value)))
            else:
                value = six.text_type(value)
                if " " in value:
                    value = "'" + value + "'"
                keywords.append("{key:s}={value:s}".format(key=key, value=value))
        return keywords
        
    def pos_string(self):
        """Position string, in a nicely formatted way."""
        try:
            position = self.position.transform_to('icrs')
        except (AttributeError, ValueError):
            string = self.position.to_string()
        else:
            string = "{ra:s} {dec:s}".format(
                ra = position.ra.to_string(u.hourangle, precision=3, pad=True),
                dec = position.dec.to_string(u.deg, precision=3, pad=True, alwayssign=True),
                )
        return string
        
    def to_starlist(self):
        """Return a starlist line."""
        position = self.position.transform_to('fk5')
        name = six.text_type(self.name).strip().replace(" ","_")
        line = "{name:<15.15s} {ra:s} {dec:s} {epoch:.0f} {keywords:s}".format(
                                name = name,
                                ra = position.ra.to_string(u.hourangle, sep=' ', precision=3, pad=True),
                                dec = position.dec.to_string(u.deg, sep=" ", precision=3, pad=True, alwayssign=True),
                                epoch = position.equinox.jyear,
                                keywords = " ".join(self._repr_keywords_()))
        return line
        
    @classmethod
    def from_starlist(cls, line):
        """Parse a single line from a starlist into the Target data structure."""
        name, position, kw = parse_starlist_line(line)
        return cls(name=name, position=position, _keywords=kw)
    

class TargetList(collections.MutableSequence, list):
    """A target list"""
    def __init__(self, iterable=None):
        super(TargetList, self).__init__()
        if iterable is not None:
            self.extend(iterable)
    
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
        return list.__setitem__(self, key, self._type_check(value))
    
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
        r = list.__getitem__(self, key)
        if isinstance(r, list):
            return self.__class__(r)
        return r
        
    def __add__(self, item):
        """Add two TargetList objects together."""
        return self.__class__(super(TargetList, self).__add__(item))
        
    def __mul__(self, value):
        """Multiply"""
        return self.__class__(super(TargetList, self).__mul__(value))
    
    def __rmul__(self, value):
        """Reverse multiply."""
        return self.__class__(super(TargetList, self).__rmul__(value))
        
    def __len__(self):
        """Length, from the underlying list."""
        return list.__len__(self)
        
    def insert(self, index, item):
        """Insert an item, and check type."""
        return list.insert(self, index, self._type_check(item))
    
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
        
    def to_starlist_stream(self, file):
        """Write to a stream"""
        for t in self:
            file.write(t.to_starlist())
            file.write("\n")
        
    def to_starlist(self, filename, mode='w'):
        """Write to a starlist file."""
        if isinstance(filename, io.IOBase) or hasattr(filename, 'write'):
            stream = filename
            self.to_starlist_stream(stream)
        else:
            with open(filename, mode) as stream:
                self.to_starlist_stream(stream)
        