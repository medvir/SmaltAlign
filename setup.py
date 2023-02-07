#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    use_scm_version=True,
    setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
    install_requires=['setuptools_scm'],
    tests_require=['pytest', 'flake8'],
    name='SmaltAlign',
    description='Calling consensus sequence',
    url='https://github.com/medvir/SmaltAlign',
    author='Maryam Zaheri, Stefan Schmutz, Michael Huber',
    author_email='zaheri@gmail.com, stefan.schmutz@uzh.ch, huber.michael@virology.uzh.ch',
    license='MIT',
    packages=find_packages('src'),  # include all packages under src
    package_dir={'': 'src'},  # tell setuptools packages are under src
    #package_data={'smaltalign': ['references/*']},
    entry_points={
        'console_scripts': ['smaltalign = smaltalign.cli:main'],
        },
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',

        # Pick your license as you wish (should match "license" above)
        'License :: MIT',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        #'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.9']
)
