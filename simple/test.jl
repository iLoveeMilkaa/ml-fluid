

using WaterLily
function circle(n,m;Re=100,U=1)
    # signed distance function to circle
    radius, center = m/8, m/2-1
    sdf(x,t) = √sum(abs2, x .- center) - radius

    Simulation((n,m),   # domain size
               (U,0),   # domain velocity (& velocity scale)
               2radius; # length scale
               ν=U*2radius/Re,     # fluid viscosity
               body=AutoBody(sdf)) # geometry
end

circ = circle(3*2^5,2^6)
sim_step!(circ)

using Plots
u = circ.flow.u[:,:,1] # first component is x
plot = contourf(u') # transpose the array for the plot
