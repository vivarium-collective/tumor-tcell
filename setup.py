import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='tumor-tcell',
    version='0.0.7',
    packages=[
        'tumor_tcell',
        'tumor_tcell.composites',
        'tumor_tcell.processes',
        'tumor_tcell.experiments',
        'tumor_tcell.library',
    ],
    author='John Hickey, Eran Agmon',
    author_email='jwhickey@stanford.edu, eagmon@stanford.edu',
    url='https://github.com/vivarium-collective/tumor-tcell',
    license='MIT',
    entry_points={
        'console_scripts': []},
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    include_package_data=True,
    install_requires=[
        'vivarium-core',
        'vivarium-cell',
        'pymunk'])
