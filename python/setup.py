#!/usr/bin/env python

import os
import glob
import platform
import re
import shutil
import subprocess
import sys
import sysconfig
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.command.clean import clean
from setuptools import setup, Extension
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.errors import DistutilsExecError
from distutils.version import StrictVersion
from distutils.util import convert_path
from distutils.sysconfig import customize_compiler
from distutils.ccompiler import show_compilers


# Path to the directory that contains this setup.py file.
base_dir = os.path.abspath(os.path.dirname(__file__))


def find_package_version(package_name):
    try:
        return __import__(package_name).__version__
    except ImportError:
        return None


def is_package_installed(package_name):
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False


def in_conda():
    return os.path.exists(os.path.join(sys.prefix, 'conda-meta'))


def find_version():
    """Extract the current version of Cryo2Grid.
    
    The version will be obtained from (in priority order):
    1. version.py (file created only when using GitHub Actions)
    2. git describe
    3. __init__.py (as failback)

    """
    version_file = os.path.join(base_dir, 'cryo2grid', 'version.py')
    if os.path.isfile(version_file):
        with open(version_file) as f:
            version = f.read().strip()

        print('Version found: %s (from version.py)' % version)
        return version

    try:
        git_output = subprocess.check_output(['git', 'describe', '--abbrev=7', '--dirty', '--always', '--tags'])
        git_output = git_output.strip().decode()

        if git_output.startswith('v'):
            git_output = git_output[1:]
        version = git_output.replace('dirty', 'mod').replace('-', '+', 1).replace('-', '.')

        print('Version found %s (from git describe)' % version)
        return version
    except:
        pass
    
    init_file = os.path.join(base_dir, 'cryo2grid', '__init__.py') 
    with open(init_file) as f:
        for line in f:
            version_match = re.match(r'^__version__ = "(.+?)"$', line)

            if version_match:
                version = version_match.group(1)

                print('Version found %s (from __init__.py)' % version)
                return version
    
    raise RuntimeError('Could not find version string for Cryo2Grid.')


def execute_command(cmd_line):
    """Simple function to execute bash command."""
    args = cmd_line.split()
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    output, errors = p.communicate()
    return output, errors


class CustomBuild(build):
    """Ensure build_ext runs first in build command."""
    def run(self):
        # Fix to make it compatible with wheel
        # We copy src directory only when it is not present already.
        # If it is present, it means that we are creating linux wheels using the manylinux docker image
        # The src copy is done outside setup.py before we start creating the wheels
        # Source: https://github.com/pypa/pip/issues/3500
        if not os.path.exists('src'):
            shutil.copytree('../include/', 'src',ignore=shutil.ignore_patterns('*.o', '*.d'))

        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    """Ensure build_ext runs first in install command."""
    def run(self):
        # This is not called when creating wheels for linux in the docker image
        if not os.path.exists('src'):
            shutil.copytree('../include', 'src',ignore=shutil.ignore_patterns('*.o', '*.d'))

        self.run_command('build_ext')
        install.run(self)

        # It means that we are in python, so we can safely remove src
        if not os.path.exists('python'):
            shutil.rmtree('src')


class CustomSdist(sdist):
    """Add swig interface files into distribution from parent directory."""
    def make_release_tree(self, base_dir, files):
        sdist.make_release_tree(self, base_dir, files)
        link = 'hard' if hasattr(os, 'link') else None
        pkg_dir = os.path.join(base_dir, 'cryo2grid')
        self.copy_file(os.path.join('cryo2grid', 'cryo2grid.i'), pkg_dir, link=link)

    def run(self):
        if not os.path.exists('src'):
            shutil.copytree('../include/', 'src',ignore=shutil.ignore_patterns('*.o', '*.d'))

        sdist.run(self)

        # It means that we are in python, so we can safely remove src
        if not os.path.exists('python'):
            shutil.rmtree('src')


class CustomBuildExt(build_ext):
    """Custom build_ext to set SWIG options and print a better error message."""
    def finalize_options(self):
        # Setting include_dirs, library_dirs, swig_opts here instead of in Extension constructor allows them to be
        # overridden using -I and -L command line options to python setup.py build_ext.
        build_ext.finalize_options(self)

        # Cryo2Grid
        self.include_dirs.append('src')
        # SWIG
        # shadow, creates a pythonic wrapper around cryo2grid
        # castmode
        self.swig_opts = ['-c++', '-small', '-naturalvar', '-fastdispatch', '-shadow']
        self.swig_opts += ['-I%s' % i for i in self.include_dirs]

        print('- include_dirs: %s\n- library_dirs: %s' % (self.include_dirs, self.library_dirs))
        print('- swig options: %s' % self.swig_opts)
        print('- libraries: %s' % self.libraries)

    def swig_sources(self, sources, extension):
        try:
            return build_ext.swig_sources(self, sources, extension)
        except DistutilsExecError:
            print('\nError: SWIG failed.',
                  'You may need to manually specify the location of include and library directories. '
                  'For example:',
                  '  python setup.py build_ext -I{} -L{}'.format(self.include_dirs, self.library_dir),
                  '  python setup.py install',
                  sep='\n')
            sys.exit(1)

    def build_extensions(self):
        customize_compiler(self.compiler)

        include_system = set(self.include_dirs).intersection(['/usr/local/include', '/usr/include'])
        # Use g++ by default
        self.compiler.compiler_so[0] = "g++"
        if platform.system() == 'Darwin':
            # To get the right @rpath on macos for libraries
            self.extensions[0].extra_link_args.append('-Wl,-rpath,' + '/usr/lib')
            # Use clang++ on macOS (use <brew install llvm> to install to get OpenMP)
            self.compiler.compiler_so[0] = "clang++"

        print('- extra link args: %s' % self.extensions[0].extra_link_args)

        remove_flags = ["-Wstrict-prototypes", "-Wall", "-Wsign-compare", "-g"]
        for remove_flag in remove_flags:
            try:
                self.compiler.compiler_so.remove(remove_flag)
            except ValueError:
#                print('Warning: compiler flag %s is not present, cannot remove it.' % remove_flag)
                pass

        self.compiler.compiler_so.insert(2, "-Wno-deprecated")
        self.compiler.compiler_so.append("-std=gnu++11")
        self.compiler.compiler_so.append("-Wno-long-long")
        self.compiler.compiler_so.append("-fopenmp")

        print('- compiler options: %s' % self.compiler.compiler_so)
        build_ext.build_extensions(self)


# Dirty fix when the compilation is not done in ./python
# Fix for readthedocs
c2g_swig_interface = 'cryo2grid/c2g.i'
package_dir = {}
if os.path.exists('python'):
    c2g_swig_interface = 'python/' + c2g_swig_interface
    package_dir = {'cryo2grid': 'python/cryo2grid'}

c2g_extension = Extension(
    'cryo2grid._c2g_wrapper',
    sources=[c2g_swig_interface],
    swig_opts=['-threads']
)

setup(
    name='cryo2grid',
    version=find_version(),
    author='Andreas F. Tillack, Althea A. Hansel, Matthew Holcomb, Stefano Forli',
    author_email='forli@scripps.edu',
    license='LGPL-2.1',
    url='https://forlilab.org',
    description='Python interface to Cryo2Grid',
    long_description=open(os.path.join(base_dir, 'README.md')).read(),
    long_description_content_type="text/markdown",
    zip_safe=False,
    cmdclass={'build': CustomBuild, 'build_ext': CustomBuildExt, 'install': CustomInstall, 'sdist': CustomSdist},
    packages=['cryo2grid'],
    package_dir=package_dir,
    python_requires='>=3.5.*',
    ext_modules=[c2g_extension],
    classifiers=[
        'Environment :: Console',
        'Environment :: Other Environment',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2.1 (LGPL-2.1)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        #'Operating System :: Microsoft :: Windows',
        'Operating System :: OS Independent',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries'
    ]
)
