from os import path
from setuptools import setup, find_packages


setup(name = 'pyCM',
	version = '0.4.0.dev0',
	description = 'Python contour method toolkit',
	long_description = 'https://github.com/majroy/pyCM',
	url = 'https://github.com/majroy/pyCM',
	author = 'M J Roy',
	author_email = 'matthew.roy@manchester.ac.uk',

	classifiers=[
		'Development Status :: 2 - Pre-Alpha',
		'Environment :: Win32 (MS Windows)',
		'Topic :: Scientific/Engineering :: Visualization',
		'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
		'Operating System :: Microsoft :: Windows :: Windows 7',
		'Programming Language :: Python :: 2.7',
		'Framework :: Spyder',
		'Intended Audience :: End Users/Desktop',
		'Natural Language :: English',
		],

	install_requires=['vtk>=6.0','numpy','scipy','pyyaml','matplotlib'],
	license = 'Creative Commons Attribution-Noncommercial-Share Alike license',
	keywords = 'residual stress contour method VTK',
	packages=['pyCM', 'pyCM.meta'],
	package_data = {'pyCM' : ['README.MD',], 'pyCM.meta' : ['pyCM_logo.png',] },
	include_package_data=True

	)
