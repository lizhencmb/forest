#!/usr/bin/env python
from setuptools import setup

setup(
    name='forest',
    version='0.1',
    py_modules=['forest'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        forest=forest:cli
    '''
)