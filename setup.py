from os import path
import sys
from setuptools import setup, find_packages


setup(name = 'pyCM',
    version = '2.0.0post1',
    description = 'Python contour method toolkit',
    long_description = 'https://github.com/majroy/pyCM',
    url = 'https://github.com/majroy/pyCM',
    author = 'M J Roy',
    author_email = 'matthew.roy@manchester.ac.uk',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Win32 (MS Windows)',
        'Topic :: Scientific/Engineering :: Visualization',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: Microsoft :: Windows :: Windows 7',
        'Programming Language :: Python :: 3.7',
        'Intended Audience :: End Users/Desktop',
        'Natural Language :: English',
        ],

    install_requires=['vtk>=6.0','numpy','scipy','pyyaml>=5.0','matplotlib','PyQt5>=5','h5py','sklearn','shapely','pyclipper'],
    license = 'Creative Commons Attribution-Noncommercial-Share Alike license',
    keywords = 'residual stress contour method VTK',
    packages=['pyCM', 'pyCM.meta'],
    package_data = {'pyCM' : ['README.MD',], 'pyCM.meta' : ['*.*',] },
    include_package_data=True

    )
