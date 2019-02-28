#!/usr/bin/env python
import sys
from setuptools import setup


if sys.version_info < (3, 5):
    raise RuntimeError('Python >= 3.5 is required')

setup()
