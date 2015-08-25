UKIRT Faint Standards Starlists
===============================

This example shows how to use KOPy to generate a starlist of the UKIRT Faint Standards, which are useful infrared faint standard stars.

Usage
-----

To make a starlist of UKIRT Faint Standards, run::
    
    $ ukirtfs_starlist.py
    

The starlist will be written to standard output. To save it to a file directly, you can use the ``-o`` option::
    
    $ ukirtfs_starlist.py -o UKIRTFS.txt
    

To make a starlist only containing the primary UKIRT faint standards, use the ``--fs`` flag::
    
    $ ukirtfs_starlist.py --fs
    

