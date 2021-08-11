
import pkg_resources, sys
from setuptools import setup
from setuptools.extension import Extension

def is_installed(requirement):
    """Checks if package is installed
    
    The argument should be given in the following form: 
    package >= version
    
    """
    
    try:
        pkg_resources.require(requirement)
    except pkg_resources.ResolutionError:
        return False
    else:
        return True

def checkDep(package, version):
    """Checks if package is installed with a minimum version
    
    If the package does not meet the requirement, a message is
    given with instructions on how to update it.
    Then, the installation is halted.
    
    """
    
    if not is_installed( package + " >= " + version ):
        import textwrap
        print(textwrap.dedent("""
            {0} >= {1} is required, but {2} is installed.
            You can upgrade {0} with:
            $ pip install -U {0}
            """.format(package, version, pkg_resources.get_distribution(package).version) ),
            file=sys.stderr)
        exit(1)
    else:
        print("Meets " + package + " >= " + version)

# Check for minimum version of setuptools    
checkDep("setuptools", "18.0")

# Check for minimum version of pip
checkDep("pip", "19.0")

# This allows us to depend on cython while installing a module
# in an environment that doesn't necessarily have cython installed yet.
try:
    from Cython.Build import build_ext
    print("Using Cython.Build.build_ext!")
except:
    # If we couldn't import Cython, use the normal setuptools
    # and look for a pre-compiled .c file instead of a .pyx file
    print("Using setuptools.command.build_ext!")
    from setuptools.command.build_ext import build_ext
    ext_modules = [Extension("odecell.paropt", ["odecell/paropt.c"])]
else:
    # If we successfully imported Cython, look for a .pyx file
    ext_modules = [Extension("odecell.paropt", ["odecell/paropt.pyx"])]

#from setuptools.command.build_ext import build_ext
#ext_modules = [ Extension("odecell.paropt", ["odecell/paropt.pyx"]) ]

# This allows cython code to depend on NumPy *and* be installed in a system
# that doesn't have NumPy yet. Setuptools will only call this class after 
# installing all dependencies, so it is safe to `import numpy` here.
class CustomBuildExtCommand(build_ext):
    """build_ext command for use when numpy headers are needed."""
    def run(self):

        # Import numpy here, only when headers are needed
        import numpy

        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())

        # Call original build_ext command
        build_ext.run(self)

def readme():
    with open('README.rst') as f:
        return f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='odecell',
    version='1.0.0',
    author='Marcelo C. R. Melo',
    author_email='melomcr@gmail.com',
    license='GPLv3',
    packages=['odecell'],
    description='Cell-Scale ODE Environment',
    long_description=readme(),
    long_description_content_type="text/markdown",
    keywords='chemistry cell-biology ODE',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3) ',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.6.8',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Pharmacokinetic',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    
    #setup_requires=requirements,
    install_requires=requirements,
    
    cmdclass = {'build_ext': CustomBuildExtCommand},
    ext_modules = ext_modules,
    
    zip_safe=False
) 

### Release

# Before making a release, run 'python setup.py build_ext' first, to ensure that paropt.c is present and up-to-date for the source code distribution.

### DEBUG

#For this purpose, the DISTUTILS_DEBUG environment variable can be set to anything except an empty string, and distutils will now print detailed information about what it is doing, dump the full traceback when an exception occurs, and print the whole command line when an external program (like a C compiler) fails.

###


