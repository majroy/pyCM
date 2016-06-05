from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.MD')) as f:
	long_description = f.read()

setup(name = 'pyCM',
	version = '0.1.0.dev0',
	description = 'Python contour method toolkit',
	long_description = long_description,
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

	install_requires=['vtk','numpy','scipy'],
	license = 'Creative Commons Attribution-Noncommercial-Share Alike license',
	keywords = 'residual stress contour method VTK',
	packages=['pyCM'],
	package_data = {'pyCM' : ['meta/*',] },

	)
