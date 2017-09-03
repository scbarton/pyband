# pyband
This is a python wrapper for BAND.F by John S. Newman with portions by others. It uses `f2py` from `numpy` and is therefore restricted to Pythonn 2.7.

Portions are copywrite John S. Newman 1998, and there are other contributions as listed in License.txt.

## To install

    $ pip install git+git://github.com/scbarton/pyband.git

## To test
Using [Jupyter Notebook](http://jupyter.org)

    %pylab inline
    from pyband.example import samplefd
    samplefd();

which should produce the following plot:

![plot](/path/to/img.jpg "Title")

You may also use the attached notebook [pyband.ipynb]() which includes a performance comparison with `scipy.integrate.solve_bvp`.
  
