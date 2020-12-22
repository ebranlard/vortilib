[![Build Status](https://travis-ci.com/ebranlard/vortilib.svg?branch=master)](https://travis-ci.com/ebranlard/vortilib)
<a href="https://www.buymeacoffee.com/hTpOQGl" rel="nofollow"><img alt="Donate just a small amount, buy me a coffee!" src="https://warehouse-camo.cmh1.psfhosted.org/1c939ba1227996b87bb03cf029c14821eab9ad91/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f446f6e6174652d4275792532306d6525323061253230636f666665652d79656c6c6f77677265656e2e737667"></a>

# vortilib

Vortex-methods and vortex-theory library, with ressources in: python, matlab, fortran, C and mathematica.


# Installation and testing
```bash
git clone http://github.com/ebranlard/vortilib
cd vortilib
python -m pip install -r requirements.txt
python -m pip install -e .
pytest
```

## Sub-packages

The vortilib package contains the following packages:
- elements: simple vortex/source/doublet elements (points, segments, cylinders, rings) functions to compute: induced velocities, streamfunction, velocity potential.  
- particles: tools for vortex particles method (interpolation, projection)
- mesh: tools to handle meshes


## References and how to cite
If you find some of this repository useful and use it in your research, thank you for using the following citation: 
```bibtex
@book{Branlard:book,
    author = {E. Branlard},
    title = {Wind Turbine Aerodynamics and Vorticity-Based Methods: Fundamentals and Recent Applications},
    year = {2017},
    publisher= {Springer International Publishing},
    doi={10.1007/978-3-319-55164-7},
    isbn={ 978-3-319-55163-0}
}
```


## Contributing
Any contributions to this project are welcome! If you find this project useful, you can also buy me a coffee (donate a small amount) with the link below:


<a href="https://www.buymeacoffee.com/hTpOQGl" rel="nofollow"><img alt="Donate just a small amount, buy me a coffee!" src="https://warehouse-camo.cmh1.psfhosted.org/1c939ba1227996b87bb03cf029c14821eab9ad91/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f446f6e6174652d4275792532306d6525323061253230636f666665652d79656c6c6f77677265656e2e737667"></a>
