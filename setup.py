from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

###### calculate.pyx ########
setup(
	cmdclass = {'build_ext': build_ext},
	ext_modules = [Extension("calculate", ["calculate.pyx"])]
)