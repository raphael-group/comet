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

from distutils.core import setup, Extension

module1 = Extension('comet',
        include_dirs=[numpy.get_include()],
        libraries = ['glib-2.0'],
        extra_compile_args = compile_args,
        sources = ['utils/cephes/polevl.c','utils/cephes/gamma.c', 'utils/cephes/incbet.c', 'weights.c', 'utils/utilities.c',
                   'mutation_data.c', 'cometmodule.c', 'comet_mcmc.c', 'comet_exhaustive.c'],
        )
setup (name = 'CoMEt',
        version = '-1',
        description = 'Runs CoMEt.',
        ext_modules = [module1]
        )
