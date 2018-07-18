from setuptools import setup
import os.path

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

pth = os.path.dirname(os.path.abspath(__file__))+ '/requirements.txt'

REQUIREMENTS = [i.strip() for i in open(pth).readlines()]

setup(name='spaghetti', #name of package
      version='1.0.1',
      description='Spatial Network Analysis (SPAtial GrapHs: nETworks, Topology, & Inference)', #short <80chr description
      url='https://github.com/pysal/spaghetti', #github repo
      maintainer='James Gaboardi',
      maintainer_email='jgaboardi@gmail.com',
      test_suite = 'nose.collector',
      tests_require=['nose'],
      keywords='spatial statistics, networks, graphs',
      classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: GIS',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6'
        ],
      license='3-Clause BSD',
      packages=['spaghetti'],
      install_requires=REQUIREMENTS,
      zip_safe=False,
      cmdclass = {'build.py':build_py})
