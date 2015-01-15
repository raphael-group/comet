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
    try:
        return out
    except:
        print "Error: " + err

glib_flgs = subprocessOutput([ 'pkg-config', '--cflags', 'glib-2.0' ]).split()
glib_libs = subprocessOutput([ 'pkg-config', '--libs', 'glib-2.0' ]).split()
compile_args = glib_flgs
compile_args += ['-g', '-O0']

srcs = ['utils/cephes/polevl.c','utils/cephes/gamma.c', 'utils/cephes/incbet.c',
        'weights.c', 'utils/utilities.c', 'mutation_data.c', 'cometmodule.c',
        'comet_mcmc.c', 'comet_exhaustive.c']
module = Extension('cComet', include_dirs=[numpy.get_include()], sources = srcs,
                   libraries = ['glib-20'], extra_compile_args = compile_args)
setup(name='CoMEt', version='1.0', description='C module for running CoMEt.',
      ext_modules=[module])

############################################################################
# Second compile the Fortran code

# Load required modules
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

# Compile the bipartite edge swap code
config = Configuration('', '', '')
config.add_extension('bipartite_edge_swap', sources=['bipartite_edge_swap.f95'])
setup(**config.todict())