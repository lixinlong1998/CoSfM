# file: setup.py
from distutils.core import setup
from Cython.Build import cythonize

setup(name='CTPsGenerator',
      ext_modules=cythonize(r"D:\Research\20221223_CoSfM\Release\CFTM_v1.2\src_CFTM\Func_CTPsGenerator.py",
                            language="3"))
