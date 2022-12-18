# poissonSolver

Numerical solution for second order nonlinear differential equation (Poisson) 

Solver expects 2 layers to be defined before start:

1. Define global simulation parameters of bias and temperature by:

*setSimParameters(bias=0.0, temp=100);*

2. First initilize a layer with **newLayer()**

3. Change layer property with *layerProp(args:

                              layerLength = 1e-6,

                              dopingDensity = 1e16,

                              dopingType = 'n',

                              materialType = 'mct',

                              materialComp = '0.3')*
