#!/usr/bin/env python

from distutils.core import setup

VERSION='1.0dev'

setup(name='consScore',
      version=VERSION,
      author='Suliat Yakubu',
      author_email='suliat16@gmail.com',
      license='MIT',
      description='A set of methods combining programs to calculate the conservation score ' 
                  'of the amino acids in a sequence',
      url= 'https://github.com/suliat16/consScore',
      packages=['consScore'],
      package_data={'example_data': []},
      classifiers=['License :: OSI Approved :: MIT License',
                   'Programming Language :: Python 3'
                   'Topic :: Scientific/Engineering :: Bio-Informatics'
                   'Operating System :: POSIX :: Linux'
                   'Intended Audience :: Science/Research'
                   'Development Status :: 2 - Pre-Alpha']
      )
