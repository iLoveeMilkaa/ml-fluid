# Solve the Navier-Stokes equations for a 2D incompressible fluid using a simple finite difference method.
import LinearAlgebra as LA
import Plots as PLT

# Gleichung:
# \partial_t v + (v \cdot \nabla) v = 1/ rho * (-\nabla p + \nu \Delta v)
