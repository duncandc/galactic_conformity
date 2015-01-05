#Duncan Campbell
#Yale University
#November, 2014
#setup build_perfect_groups.pyx

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

ext_modules = [Extension("host_props", ["host_props.pyx"])]
extra_compile_args=["-O3"]

setup(
  name = 'host_props app',
  cmdclass = {'build_ext': build_ext},
  include_dirs=[numpy.get_include()],
  ext_modules = ext_modules
)

#to compile code type:
#    python setup_host_props.py build_ext --inplace