
from setuptools import setup
from distutils.command.build_py import build_py


package = 'spaghetti'

# Get __version__ from package/__init__.py
with open(package+'/__init__.py', 'r') as f:
    exec(f.readline())

description = 'Analysis of Network-constrained Spatial Data'

# Fetch README.md for the `long_description`
with open('README.md', 'r', encoding='utf-8') as file:
    long_description = file.read()


def _get_requirements_from_files(groups_files):
    """returns a dictionary of all requirements
    keyed by type of requirement.
    
    Parameters
    ----------
    
    groups_files : dict
        k - descriptive name, v - file name (including extension)
    
    Returns
    -------
    
    groups_reqlist : dict
        k - descriptive name, v - list of required packages
    """
    groups_reqlist = {}
    for k,v in groups_files.items():
        with open(v, 'r') as f:
            pkg_list = f.read().splitlines()
        groups_reqlist[k] = pkg_list
    return groups_reqlist


def setup_package():
    """sets up the python package"""
    
    # Requirements for: base, dev, docs, plus, and test builds
    _groups_files = {'base': 'requirements.txt',
                     'dev': 'requirements_dev.txt',
                     'docs': 'requirements_docs.txt',
                     'plus': 'requirements_plus.txt',
                     'tests': 'requirements_tests.txt'}
    reqs = _get_requirements_from_files(_groups_files)
    install_reqs = reqs.pop('base')
    extras_reqs = reqs

    setup(name=package,
          version=__version__,
          description=description,
          long_description = long_description,
          long_description_content_type='text/markdown',
          url='https://github.com/pysal/'+package,
          download_url='https://pypi.org/project/'+package,
          maintainer='James D. Gaboardi',
          maintainer_email='jgaboardi@gmail.com',
          test_suite = 'nose.collector',
          tests_require=['nose'],
          keywords='spatial statistics, networks, graphs',
          classifiers=['Development Status :: 5 - Production/Stable',
                       'Intended Audience :: Science/Research',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Education',
                       'Topic :: Scientific/Engineering',
                       'Topic :: Scientific/Engineering :: GIS',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python',
                       'Programming Language :: Python :: 3.6',
                       'Programming Language :: Python :: 3.7'],
          license='3-Clause BSD',
          packages=[package],
          py_modules=[package],
          install_requires=install_reqs,
          extras_require=extras_reqs,
          zip_safe=False,
          cmdclass = {'build.py':build_py},
          python_requires='>3.5')


if __name__ == '__main__':
    
    setup_package()
