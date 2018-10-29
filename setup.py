#!/usr/bin/env python

import setuptools

with open('README.md', 'r') as fh:
      LONG_DESCRIPTION = fh.read()
VERSION='1.0.0'

setuptools.setup(name='consScore',
      version=VERSION,
      author='Suliat Yakubu',
      author_email='suliat16@gmail.com',
      license='MIT',
      description='A set of methods combining programs to calculate the conservation score ' 
                  'of the amino acids in a sequence',
      long_description = LONG_DESCRIPTION,
      long_description_content_type="text/markdown",
      url= 'https://github.com/suliat16/consScore',
      packages= setuptools.find_packages(),
      include_package_data=True,
      classifiers=['License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Operating System :: POSIX :: Linux',
                   'Intended Audience :: Science/Research',
                   'Development Status :: 2 - Pre-Alpha'],
      )