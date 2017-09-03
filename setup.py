from numpy.distutils.core import setup, Extension
import setuptools

ext = Extension(name = 'band', sources = ['BAND.pyf', 'BAND.F'])
                 
if __name__ == "__main__":
    setup(name = 'pyband',
          version='1.0.0a1',
          description       = "python wrapper for BAND.F",
          author            = "Scott Calabrese Barton",
          author_email      = "scb@msu.edu",
          ext_package='pyband',
          ext_modules = [ext],
          install_requires=['numpy'],
          python_requires='~=2.7',
          py_modules=['pyband'],
          packages=setuptools.find_packages(),
          )

