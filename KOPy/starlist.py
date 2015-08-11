"""
This is a parser for the Keck starlist format.
The format is documented at the Keck website: <https://www2.keck.hawaii.edu/observing/starlist.html>

Author: Alex Rudy, UCSC
Date: 2015-01-24
License: BSD
"""

import warnings
from datetime import date, datetime
import io
import contextlib
import astropy.units as u
import astropy.time
from astropy.coordinates import SkyCoord, FK4, FK5, AltAz
from astropy.utils.data import get_readable_fileobj
from collections import OrderedDict

from . import version

import re
_starlist_re_raw = r"""
    ^(?P<Name>.{1,15})[\s]+ # Target name must be the first 15 characters.
    (?P<RA>(?:[\d]{1,2}[\s:h][\s\d]?[\d][\s:m][\s\d]?[\d](?:\.[\d]+)?s?)|(?:[\d]+\.[\d]+))[\s]+  # Right Ascension, HH MM SS.SS+
    (?P<Dec>(?:[+-]?[\d]{1,2}[\s:d][\s\d]?[\d][\s:m][\s\d]?[\d](?:\.[\d]+)?s?)|(?:[\d]+\.[\d]+)) # Declination, (-)DD MM SS.SS+
    (?:[\s]+(?P<Equinox>(?:[\d]{4}(?:\.[\d]+)?)|(?:[A-Za-z]+)))?[\s]* # Equinox.
    (?P<Keywords>.+)?$ # Everything else must be a keyword.
    """
    
_starlist_re = re.compile(_starlist_re_raw, re.VERBOSE)

_starlist_re_strict = r"""
    ^(?P<Name>.{15})\  # Target name must be the first 15 characters.
    (?P<RA>[\d]{2}\ [\d]{2}\ [\d]{2}(?:\.[\d]+)?)\   # Right Ascension, HH MM SS.SS+
    (?P<Dec>[+-]?[\d]{1,2}\ [\d]{2}\ [\d]{2}(?:\.[\d]+)?)\  # Declination, (-)DD MM SS.SS+
    (?P<Equinox>(?:[\d]{4}(?:\.[\d]+)?)|(?:APP))[\ ]? # Equinox.
    (?P<Keywords>[^\ ].+)?$ # Everything else must be a keyword.
    """

_starlist_token_parts = ["Name", "RA", "Dec", "Equinox", "Keywords"]

def tokenize(text, identifier="<stream>"):
    """Parse a starlist line into individual tokens, looking for the first bad token."""
    token_res = [re.compile(t) for t in _starlist_re_raw.splitlines()]
    for token_re, token_name in zip(token_res, _starlist_token_parts):
        m = token_re.match(text)
        if not m:
            raise ValueError("Couldn't parse token {:s} from '{:s}'".format(token_name, text))
        text = text[m.start():m.end()]
        

def verify_starlist_line(text, identifier="<stream>", warning=False):
    """Verify that the given line is a valid starlist."""
    line = text
    messages = []
    match = _starlist_re.match(text)
    if not match:
        raise ValueError("Couldn't parse '{0:s}', no regular expression match found.".format(text))
    data = match.groupdict("")
    
    # Check the Name:
    name_length = match.end('Name') - match.start('Name') + 1
    if name_length < 15:
        messages.append(('Warning','Name','Name should be exactly 15 characters long (whitespace is ok.) len(Name)={0:d}'.format(name_length)))
    
    # Check the RA starting position.
    if match.start('RA') + 1 != 17:
        messages.append(('Error','RA','RA must start in column 17. Start: {0:d}'.format(match.start('RA')+1)))
    
    # Check the Dec starting token
    if match.start('Dec') - match.end('RA') != 1:
        messages.append(('Warning','Dec','RA and Dec should be separated by only a single space, found {0:d} characters.'.format(match.start('Dec') - match.end('RA'))))
    
    # Check the Equinox starting token.
    if match.start('Equinox') - match.end('Dec') != 1:
        messages.append(('Warning','Equinox','Dec and Equinox should be separated by only a single space, found {0:d} characters.'.format(match.start('Equinox') - match.end('Dec'))))
    
    if match.group("Keywords") and len(match.group("Keywords")):
        for kwarg in match.group("Keywords").split():
            if kwarg.count("=") < 1:
                messages.append(('Error', 'Keywords', 'Each keyword/value pair must have 1 "=", none found {!r}'.format(kwarg)))
            if kwarg.count("=") > 1:
                messages.append(('Error', 'Keywords', 'Each keyword/value pair must have 1 "=", {0:d} found {1!r}'.format(kwarg.count("="), kwarg)))
    
    composed_messages = []
    for severity, token, message in messages:
        composed_message = "[{0}][{1} {2}] {3}".format(severity, identifier, token, message)
        if warning:
            warnings.warn(composed_message)
        else:
            composed_messages.append(composed_message)
    return composed_messages

def parse_starlist_line(text):
    """Parse a single line from a Keck formatted starlist, returning a dictionary of parsed values.
    
    :param text: The starlist text line.
    :returns: A dictionary of starlist object properties, set from teh starlist line.
    :raises: ValueError if the line couldn't be parsed.
    
    This function parses a single line from a starlist and returns a dictionary of items from that line. The followig keys are included:
    - `Name`: The target name.
    - `Position`: An astropy.coordinates object representing this position.
    - Any other keyword/value pairs, which are found at the end of the starlist line, and formatted as ``keyword=value``
    
    """
    match = _starlist_re.match(text)
    if not match:
        raise ValueError("Couldn't parse '{}', no regular expression match found.".format(text))
    data = match.groupdict("")
    if data.get('Equinox','') == '':
        equinox = astropy.time.Time.now()
        frame = AltAz
    elif data['Equinox'] == "APP":
        equinox = astropy.time.Time.now()
        frame = 'fk5'
    else:
        equinox = astropy.time.Time(float(data['Equinox']), format='jyear', scale='utc')
        if float(data['Equinox']) <= 1950:
            equinox = astropy.time.Time(float(data['Equinox']), format='byear', scale='utc')
            frame = 'fk4'
        else:
            frame = 'fk5'
    
    position = SkyCoord(data["RA"], data["Dec"], unit=(u.hourangle, u.degree), equinox=equinox, frame=frame)
    
    results = OrderedDict()
    for keywordvalue in data.get("Keywords","").split():
        if keywordvalue.count("=") < 1:
            warnings.warn("Illegal Keyword Argument: '{}'".format(keywordvalue))
            continue
        keyword, value = keywordvalue.split("=",1)
        keyword = keyword.strip()
        results[keyword] = value.strip()
    return data['Name'].strip(), position, results
    
def read_skip_comments(filename, comments="#"):
    """Read a filename, yielding lines that don't start with comments."""
    with get_readable_fileobj(filename) as stream:
        for line in stream_skip_comments(stream, comments):
            yield line

def stream_skip_comments(stream, comments="#"):
    """Skip comment lines from a stream."""
    for line in stream:
        if not line.startswith(comments):
            yield line.strip("\n\r").strip()
    
def parse_starlist(starlist):
    """Parse a starlist into a sequence of dictionaries."""
    for line in read_skip_comments(starlist):
        yield parse_starlist_line(line)
    
def main():
    """Command-line interface for starlist parsing."""
    import argparse
    parser = argparse.ArgumentParser(description="A Keck starlist parsing and verification tool")
    parser.add_argumnet("starlist", metavar="starlist.txt", help="starlist filename", type=argparse.FileType('r'), default='starlist.txt')
    opt = parser.parse_args()
    print("Starlist Lint {0:s} for '{1:s}'".format(version, opt.name))
    with get_readable_fileobj(starlist) as f:
        for n, line in enumerate(f):
            if not f.startswith("#"):
                messages = verify_starlist_line(line, identifier="{!r} line {:d}".format(opt.starlist, n))
                if len(messages):
                    print(line)
                    print("\n".join(messages))
    
    
