#
# setup.py
#
# Installation script to get setuptools to install VGA into
# a Python environment.
#

import sys
import setuptools

# Import the lengthy rich-text README as the package's long
# description:
with open('README.rst', 'r') as fh:
	long_description = fh.read()

setuptools_info = {
	'name': 'VGA',
	'version': '1.1.2',
	'author': 'Vlachos Research Group',
	'author_email': 'vlachos@udel.edu',
	'description': 'Python package implements the Group Additivity (GA) method for estimating thermodynamic properties',
	'long_description': long_description,
	'zip_safe': False,
	'url': 'https://github.com/VlachosGroup/VlachosGroupAdditivity',
	'packages': setuptools.find_packages(),
	'include_package_data': True,
	'exclude_package_data': {
		'': [ 'README.rst', 'docs', 'example', 'tests', 'PKG-INFO', 'LICENSE.md' ]
	    },
	'install_requires': [
		'scipy>=1.1.0',
		'numpy>=1.15.1',
		'IPython>=7.0.0',
		'PyYAML>=3.0',
	    ],
	'dependency_links': [
		'git+https://github.com/rdkit/rdkit/',
	    ],
	'classifiers': [
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering :: Chemistry",
	    ],
    }

if sys.version_info[0] >= 3:
	#
	# Augment for Python 3 setuptools:
	#
	setuptools_info['long_description_content_type'] = 'text/x-rst'

setuptools.setup(**setuptools_info)