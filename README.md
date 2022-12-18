# poissonSolver

Numerical solution for second order nonlinear differential equation (Poisson) 

Solver expects 2 layers to be defined before start:

First initilize a layer with ***newLayer()***
Change layer property with *layerProp(args:
                              layerLength = 1e-6,
                              dopingDensity = 1e16,
                              dopingType = 'n',
                              materialType = 'mct',
                              materialComp = '0.3')*
