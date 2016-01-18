#!/usr/bin/env python
"""Setup script for building leaf's python bindings"""

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(parent_package,top_path)
    config.add_extension('leaf', ['leaf/dataSpec_P5B.f90', \
        leaf/getdata.py","leaf/__init__.py"] )
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    # Global variables for this extension:
    name         = "leaf"  # name of the generated python extension (.so)
    description  = "PROSPECT-based leaf scattering method"
    long_description = "The PROSPECT RT model with flexible inputs."
    author       = "Prof. P. Lewis/NCEO & University College London"
    author_email = "p.lewis@ucl.ac.uk"
    url = "https://github.com/UCL-EO/leaf.git"
    
    setup( name=name,\
        description=description, \
        author=author, \
        author_email = author_email, \
        configuration = configuration, version="1.0.0",\
        packages=["leaf"])
