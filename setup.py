import os
from setuptools import setup

def read(fname):
	return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name = "pyCM",
	version = "1.0",
	description = "Python and VTK contour method toolkit",
	author = "M J Roy",
	author_email = "matthew.roy@manchester.ac.uk",
	url = "http://www.contour-method.manchester.ac.uk",
	keywords = "residual stress contour method VTK",
	packages=['pyCM'],
	package_data = {'pyCM' : ['meta/*',] },
	long_description = read('README'),
	classifiers=[
		"Development Status :: 2 - Pre-Alpha",
		"Environment :: Win32 (MS Windows)",
		"Topic :: Scientific/Engineering :: Visualization",
		"License :: Other/Proprietary License",
		"Operating System :: Microsoft :: Windows :: Windows 7",
		"Programming Language :: Python :: 2.7",
		"Framework :: IDLE",
		"Intended Audience :: End Users/Desktop",
		"Natural Language :: English",
		],
	)
