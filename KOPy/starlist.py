# Author: Alex Rudy, UCSC
# Date: 2015-01-24
# License: BSD
"""
This is a parser for the Keck starlist format.
The format is documented at the Keck website: https://www2.keck.hawaii.edu/observing/starlist.html

This module is a functional interface to starlist parsing. A more user-friendly object-oriented user interface is provided in :mod:`~KOPy.targets`
"""

import warnings
from datetime import date, datetime
import io
import contextlib
import six
import astropy.units as u
import astropy.time
from astropy.coordinates import SkyCoord, FK4, FK5, AltAz
from astropy.utils.data import get_readable_fileobj
from collections import OrderedDict

from . import __version__
import re

__all__ = ['tokenize', 'verify_starlist_line', 'parse_starlist_line', 'read_skip_comments', 'stream_skip_comments', 'parse_starlist', 'format_starlist_line']

_starlist_re_raw = r"""
    ^(?P<Name>.{1,15})[\s]+ # Target name must be the first 15 characters.
    (?P<RA>(?:[\d]{1,2}[\s:h][\s\d]?[\d][\s:m][\s\d]?[\d](?:\.[\d]+)?s?)|(?:[\d]+\.[\d]+))[\s]+  # Right Ascension, HH MM SS.SS+
    (?P<Dec>(?:[+-]?[\d]{1,2}[\s:d][\s\d]?[\d][\s:m][\s\d]?[\d](?:\.[\d]+)?s?)|(?:[\d]+\.[\d]+)) # Declination, (-)DD MM SS.SS+
    (?:[\s]+(?P<Equinox>(?:[\d]{4}(?:\.[\d]+)?)|(?:[A-Za-z]+)))?[\s]* # Equinox.
    (?P<Keywords>[^\#]+)? # Everything else must be a keyword.
    (?P<Comments>\#.+)?$ # Comments allowed at the end of lines.
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
    """Parse a starlist line into individual tokens, looking for the first bad token.
    
    This function will parse individual tokens with the starlist regular expression, and will raise an :exc:`ValueError` when it encounters an invalid token. 
    
    Parameters
    ----------
    text : string
        The string to tokenize.
    identifier : string
        An identifier (usually the filename) to use in the error message for when parsing fails.
    
    Raises
    ------
    ValueError
        When a token can't be parsed.
    
    """
    token_res = [re.compile(t) for t in _starlist_re_raw.splitlines()]
    for token_re, token_name in zip(token_res, _starlist_token_parts):
        m = token_re.match(text)
        if not m:
            raise ValueError("Couldn't parse token {:s} from '{:s}'".format(token_name, text))
        text = text[m.start():m.end()]
        

def verify_starlist_line(text, identifier="<stream>", warning=False):
    """Verify that the given line is a valid starlist.
    
    Parameters
    ----------
    text : string
        The string to tokenize.
    identifier : string
        An identifier (usually the filename) to use in the error message for when parsing fails.
    warning : bool
        Whether to emit the messages indicating problems with the starlist format as warnings.
    
    Raises
    ------
    ValueError
        When a line can't be parsed, even with the forgiving starlist regular expression.
        
    Returns
    -------
    messages : list
        A list of messages indicating problems with the starlist line verification.
    
    """
    line = text
    messages = []
    match = _starlist_re.match(text)
    if not match:
        raise ValueError("Couldn't parse '{0:s}', no regular expression match found.".format(text.strip("\n\r")))
    data = match.groupdict("")
    
    # Check the Name:
    name_length = match.end('Name') - match.start('Name') + 1
    if name_length < 15:
        messages.append(('WARNING','Name','Name should be exactly 15 characters long (whitespace is ok.) len(Name)={0:d}'.format(name_length)))
    
    # Check the RA starting position.
    if match.start('RA') + 1 != 17:
        messages.append(('ERROR','RA','RA must start in column 17. Start: {0:d}'.format(match.start('RA')+1)))
    #
    # # Check the Dec starting token
    # if match.start('Dec') - match.end('RA') != 1:
    #     messages.append(('WARNING','Dec','RA and Dec should be separated by only a single space, found {0:d} characters.'.format(match.start('Dec') - match.end('RA'))))
    #
    # # Check the Equinox starting token.
    # if match.start('Equinox') - match.end('Dec') != 1:
    #     messages.append(('WARNING','Equinox','Dec and Equinox should be separated by only a single space, found {0:d} characters.'.format(match.start('Equinox') - match.end('Dec'))))
    
    if match.group("Keywords") and len(match.group("Keywords")):
        for kwarg in match.group("Keywords").split():
            if kwarg.count("=") < 1:
                messages.append(('ERROR', 'Keywords', 'Each keyword/value pair must have 1 "=", none found {!r}'.format(kwarg)))
            if kwarg.count("=") > 1:
                messages.append(('ERROR', 'Keywords', 'Each keyword/value pair must have 1 "=", {0:d} found {1!r}'.format(kwarg.count("="), kwarg)))
    
    composed_messages = []
    for severity, token, message in messages:
        composed_message = "[{0:s}] {3} [{1} {2}]".format(severity, identifier, token, message)
        if warning:
            warnings.warn(composed_message)
        else:
            composed_messages.append(composed_message)
    return composed_messages

def parse_starlist_line(text):
    """Parse a single line from a Keck formatted starlist, returning a dictionary of parsed values.
    
    This uses the forgiving starlist parser, which should be robust to various errors in starlist file formats.
    
    Parameters
    ----------
    text : string
        The starlist text line.
        
    Raises
    ------
    ValueError :
        Raised if the line couldn't be parsed by the forgiving parser.
    
    Returns
    -------
    name : string
        The target name
    position : :class:`~astropy.coordinates.SkyCoord`
        The target position, as a :class:`~astropy.coordinates.SkyCoord` object.
    keywords : OrderedDict
        An ordered dictionary of keyword values applied to the starlist line.
    
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
        for expression in PARSE_KEYWORDS:
            if re.match(expression, keyword):
                try:
                    results[keyword] = PARSE_KEYWORDS[expression][1](value)
                    break
                except ValueError:
                    pass
        else:
            results[keyword] = value.strip().replace("=","")
    return data['Name'].strip(), position, results
    
def read_skip_comments(filename, comments="#"):
    """Read a filename, yielding lines that don't start with comments.
    
    Parameters
    ----------
    filename : string or filobj
        The file to be opened and read from.
    comments : string
        The string used to match the beginning of comment lines.
    
    Yields
    ------
    line : string
        Lines from the file which don't start with the comment string.
    
    """
    with get_readable_fileobj(filename) as stream:
        for line in stream_skip_comments(stream, comments):
            yield line

def stream_skip_comments(stream, comments="#"):
    """Skip comment lines from a stream.
    
    Parameters
    ----------
    stream : io.BaseIO
        The stream object to be read from.
    comments : string
        The string used to match the beginning of comment lines.
    
    Yields
    ------
    line : string
        Lines from the file which don't start with the comment string.
    
    """
    for line in stream:
        if not line.startswith(comments) and not re.match(r"^[\s]*$", line.strip("\n\r")):
            yield line.strip("\n\r").strip()
    
def parse_starlist(starlist):
    """Parse a full starlist file into a generator of target objects.
    
    Parameters
    ----------
    starlist : string or filobj
        The file to be opened and read from.
    
    Yields
    ------
    name : string
        The target name
    position : :class:`~astropy.coordinates.SkyCoord`
        The target position, as a :class:`~astropy.coordinates.SkyCoord` object.
    keywords : OrderedDict
        An ordered dictionary of keyword values applied to the starlist line.
    
    """
    for line in read_skip_comments(starlist):
        yield parse_starlist_line(line)
    
def format_starlist_position(position):
    """Output a SkyCoord object in the starlist format."""
    try:
        position = position.transform_to('fk5')
        if position.frame.equinox.jyear <= 1950:
            position = position.transform_to('fk4')
    except (ValueError, AttributeError):
        if not isinstance(position.frame, AltAz):
            raise ValueError("Can't understand frame for coordinates, should be transformable to FK5 or AltAz, got {!r}".format(position.frame))
        position_string = "{0:2.5f} {1:2.5f}".format(position.az.to(u.deg).value, position.alt.to(u.deg))
    else:
        position_string = "{ra:s} {dec:s} {epoch:.0f}".format(
            ra = position.ra.to_string(u.hourangle, sep=' ', precision=3, pad=True),
            dec = position.dec.to_string(u.deg, sep=" ", precision=3, pad=True, alwayssign=True),
            epoch = position.equinox.jyear
        )
    return position_string
    
def _format_rotator_mode(value):
    """Format rotator mode, and rais appropriate error if it can't be formatted."""
    modes = set(['pa', 'vertical', 'stationary'])
    if value.lower() not in modes:
        raise ValueError("Rotator mode must be in {!r}".format(modes))
    return value.lower()
    
    
# secondangle is one second of time in angle units.
_secondangle = u.hourangle / 3600

# Dictionary which maps regular expressions which might match keyword values to quantity transformations which 
# can output or parse values.
PARSE_KEYWORDS = {
    r'.*mag' : (lambda value : "{:.2f}".format(u.Quantity(value, u.mag).value), 
        lambda value : u.Quantity(float(value), u.mag)),
    r'pmdec': (lambda value : "{:.5f}".format(u.Quantity(value, u.arcsec / u.year).value), 
        lambda value : u.Quantity(float(value), u.arcsec/u.year)),
    r'pmra' : (lambda value : "{:.5f}".format(u.Quantity(value, _secondangle / u.year).value), 
        lambda value : u.Quantity(float(value), _secondangle / u.year)),
    r'dra': (lambda value : "{:.4}".format((u.Quantity(value,  _secondangle / u.hr)).value),
        lambda value : u.Quantity(float(value), _secondangle / u.hr)),
    r'ddec': (lambda value : "{:.3}".format(u.Quantity(value, u.arcsec/u.hr).value),
        lambda value : u.Quantity(float(value), u.arcsec/u.hr)),
    r'rotdest': (lambda value : "{:.2f}".format(u.Quantity(value, u.degree).value),
        lambda value : u.Quantity(float(value), u.degree)),
    r'rotmode': (_format_rotator_mode, lambda value : value.lower()),
    r'raoffset' : (lambda value : "{:.1f}".format(u.Quantity(value, u.arcsecond).value),
        lambda value : u.Quantity(float(value), u.arcsecond)),
    r'decoffset' : (lambda value : "{:.1f}".format(u.Quantity(value, u.arcsecond).value),
        lambda value : u.Quantity(float(value), u.arcsecond)),
}

def format_keywords(keywords):
    """Format starlist keywords in a reasonable way.
    
    This method ensures that keywords are correctly formatted for starlist output.
    It handles keywords with Quantity values and units, and ensures that all keywords
    are presented to a uniform precision if possible. When this isn't possible, the raw
    value of the keyword is converted to a string and included in the output.
    
    Parameters
    ----------
    keywords : dict-like
        A dictionary or ordered dictionary of keyword value pairs.
    
    Returns
    -------
    A list of strings, formatted for each keyword individually.
    
    """
    if not len(keywords):
        return [""]
    output = []
    output_fkeywords = []
    for key, value in keywords.items():
        for exp in PARSE_KEYWORDS:
            try:
                if re.match(exp, key):
                    output_fkeywords.append("{key}={value:s}".format(key=key, value=PARSE_KEYWORDS[exp][0](value)))
                    break
            except (ValueError, TypeError):
                pass
        else:
            value = six.text_type(value)
            if " " in value:
                value = "'" + value + "'"
            output.append("{key:s}={value:s}".format(key=key, value=value))
    return output_fkeywords + output
    
def format_starlist_line(name, position, keywords, remove_spaces=False):
    """Output the starlist."""
    name = six.text_type(name).strip()
    if remove_spaces:
        name = name.replace(" ","_")
    line = "{name:<15.15s} {position:s} {keywords:s}".format(
                            name = name,
                            position = format_starlist_position(position),
                            keywords = " ".join(format_keywords(keywords)))
    return line
    
def main():
    """Command-line interface for starlist parsing and verification."""
    import argparse
    import sys
    from astropy.utils.console import color_print
    parser = argparse.ArgumentParser(description="A Keck starlist parsing and verification tool", epilog="Parsing will be done in the 'lenient mode', with problems emitted to stderr. A correctly formatted starlist for each line, when available, will be printed to stdout, so that output can be piped into a clean starlist file.")
    parser.add_argument("starlist", metavar="starlist.txt", help="starlist filename", type=argparse.FileType('r'), default='starlist.txt')
    parser.add_argument("-o", dest='output', help="output filename", type=argparse.FileType("w"), default="-")
    opt = parser.parse_args()
    
    n_messages = 0
    all_messages = []
    for n, line in enumerate(opt.starlist):
        if line.startswith("#") or len(line.strip()) == 0:
            opt.output.write(line)
            opt.output.flush()
        else:
            identifier = "{!r} line {:d}".format(opt.starlist.name, n+1)
            try:
                messages = verify_starlist_line(line, identifier=identifier)
            except ValueError as e:
                messages = ["[ERROR] {1} [{0}]".format(identifier, e)]
            if len(messages):
                all_messages.append((n, line, messages))
            n_messages += len(messages)
            try:
                formatted_line = format_starlist_line(*parse_starlist_line(line)) + "\n"
            except ValueError:
                formatted_line = line
                opt.output.write("# WARNING Starlist lint couldn't parse next line.\n")
            opt.output.write(formatted_line)
            opt.output.flush()
            
    color_print("Starlist Lint {0:s}".format(__version__), 'green', file=sys.stderr, end="")
    sys.stderr.write(" for '{1:s}'\n".format(__version__, opt.starlist.name))
    sys.stderr.flush()
    if not len(all_messages):
        color_print("No problems found.", 'green', file=sys.stderr)
    else:
        color_print("{0:d} problems found.".format(len(all_messages)), 'yellow', file=sys.stderr)
    for n, line, messages in all_messages:
        color_print("[line {0:d}] ".format(n), 'cyan', file=sys.stderr, end="")
        color_print("=>", 'blue', file=sys.stderr, end="")
        sys.stderr.write(" '{}'\n".format(line.strip("\n")))
        for message in messages:
            if message.startswith("[ERROR]"):
                color_print(message[:len("[ERROR]")], 'red', file=sys.stderr, end="")
                sys.stderr.write(message[len("[ERROR]"):])
                sys.stderr.write("\n")
            elif message.startswith("[WARNING]"):
                color_print(message[:len("[WARNING]")], 'yellow', file=sys.stderr, end="")
                sys.stderr.write(message[len("[WARNING]"):])
                sys.stderr.write("\n")
            else:
                sys.stderr.write(message)
                sys.stderr.write("\n")
        sys.stderr.flush()
    return n_messages
