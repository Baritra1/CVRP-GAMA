from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    name="sa",
    ext_modules=cythonize(
        Extension(
            "sa",
            sources=["sa.pyx","cpu_sa.cpp"],
            include_dirs=[numpy.get_include()],
            language="c++",
            extra_compile_args=["/openmp"],
        ),
        language_level=3
    )
)