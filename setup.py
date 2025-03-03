#!/usr/bin/env python

"""
Call `pip install -e .` to install package locally for testing.
"""

from setuptools import setup

# build command
setup(
    name="miniproject",
    version="0.0.1",
    author="Selina Chen",
    author_email="yc4635@columbia.edu",
    description="A package for ...",
    classifiers=["Programming Language :: Python :: 3"],
    entry_points={
        "console_scripts": ["miniproject = miniproject.__main__:main"]
    },
)
