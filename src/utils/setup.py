'''
Created on Jul 9, 2021

@author: Lukas Graf (D-USYS, ETHZ)
'''

import setuptools
from setuptools import find_packages
from os import path

home = path.abspath(path.dirname(__file__))
with open(path.join(home, 'README.md'), encoding='utf-8') as readme:
    long_description = readme.read()


setuptools.setup(
    name='rtm_inv',
    setup_requires=['setuptools_scm'],
    use_scm_version={'version_scheme': 'post-release'},
    description = '',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    author='Lukas Graf (D-USYS, ETH Zurich)',
    author_email ='NA',
    url='NA',
    packages = find_packages(
        where = 'src'
    ),
    package_dir={'':'src'},
    include_package_data=True,
    package_data={'': []},
    classifiers = [
     "Programming Language :: Python :: 3",
     "Operating System :: OS Independent",
    ],
    install_requires=['pandas',
                      'numpy',
                      'matplotlib', 
                      'prosail',
                      'spectral'
                      ]
)
