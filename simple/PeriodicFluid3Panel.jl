import WaterLily
import Plots
import FFTW
import Random

include("PeriodicFluid.jl")   # provides random_velocity_field and make_sim

# ---------------------------------------------------------------------------
# Cell-centre velocity components and magnitude
# ---------------------------------------------------------------------------
function velocity_fields(sim)
    u = sim.flow.u
    n = size(u, 1) - 2
    uc = [(0.5 * (u[i, j, 1] + u[i+1, j, 1])) for i in 2:n+1, j in 2:n+1]
    vc = [(0.5 * (u[i, j, 2] + u[i, j+1, 2])) for i in 2:n+1, j in 2:n+1]
    mag = sqrt.(uc .^ 2 .+ vc .^ 2)
    return mag, uc, vc
end

# ---------------------------------------------------------------------------
# Run and animate — three panels side by side
# ---------------------------------------------------------------------------
sim   = make_sim(L=64, Re=1000.0)
t_end = 20.0
step  = 0.5

anim = Plots.@animate for t in range(0.0, t_end; step)
    WaterLily.sim_step!(sim, t; remeasure=false)
    mag, uc, vc = velocity_fields(sim)
    tstr = "t=$(round(WaterLily.sim_time(sim), digits=2))"

    p1 = Plots.heatmap(mag'; color=:viridis, clims=(0, 1),
                       aspect_ratio=:equal, title="|u|  $tstr",
                       xlabel="x", ylabel="y", colorbar_title="|u|")

    p2 = Plots.heatmap(uc'; color=:RdBu, clims=(-1, 1),
                       aspect_ratio=:equal, title="u_x  $tstr",
                       xlabel="x", ylabel="y", colorbar_title="u_x")

    p3 = Plots.heatmap(vc'; color=:RdBu, clims=(-1, 1),
                       aspect_ratio=:equal, title="u_y  $tstr",
                       xlabel="x", ylabel="y", colorbar_title="u_y")

    Plots.plot(p1, p2, p3; layout=(1, 3), size=(1200, 380))
end

Plots.gif(anim, "plots/periodic_fluid_3panel.gif"; fps=10)
