#!/usr/bin/python

"""Compiles the C and Fortran modules used by CoMEt."""

############################################################################
# First compile the C code

# Load required modules
from distutils.core import setup, Extension
import subprocess, numpy

def subprocessOutput(args):
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    out, err = proc.communicate()
    try: return out
    except: print "Error: " + err

glib_flgs = subprocessOutput([ 'pkg-config', '--cflags', 'glib-2.0' ]).split()
glib_libs = subprocessOutput([ 'pkg-config', '--libs', 'glib-2.0' ]).split()
compile_args = glib_flgs
compile_args += ['-g', '-O0']

srcs = ['src/c/utils/cephes/polevl.c','src/c/utils/cephes/gamma.c',
        'src/c/utils/cephes/incbet.c', 'src/c/utils/utilities.c',
        'src/c/weights.c', 'src/c/mutation_data.c', 'src/c/cometmodule.c',
        'src/c/comet_mcmc.c', 'src/c/comet_exhaustive.c']
module = Extension('cComet', include_dirs=[numpy.get_include()], sources = srcs,
                   libraries = ['glib-2.0'], extra_compile_args = compile_args)
setup(name='CoMEt', version='1.0', description='C module for running CoMEt.',
      ext_modules=[module])

############################################################################
# Second compile the Fortran code

# Load required modules
from numpy.distutils.core import Extension, setup

# Compile the bipartite edge swap code
ext = Extension(name='permute_matrix', sources=['src/fortran/permute_matrix.f95'])
setup(name='permute_matrix', ext_modules=[ext])