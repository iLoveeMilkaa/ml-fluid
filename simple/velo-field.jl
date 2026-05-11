using WaterLily
using Plots
using Dates

# Domain
const Nx, Ny = 40, 40
const T = Float32

# typsichere IC-Funktion
function make_u0(Ny::Int, ::Type{T}=Float32) where {T<:AbstractFloat}
    y0 = T(Ny) / T(2)
    w = T(Ny) / T(8)
    return (i, x) -> i == 1 ? exp(-((x[2] - y0) / w)^2) : zero(T)
end

u0 = make_u0(Ny, T)

sim = Simulation((Nx, Ny), (T(0), T(0)), Ny; U=T(1),
    ν=T(0.01),
    uλ=u0,
    perdir=(1, 2),
    body=WaterLily.NoBody(),
    T=T
)

case_id = "velo-field"
WaterLily.logger(case_id)  # schreibt Solver-Residuals nach "velo-field.log"

runtime_log = "$(case_id)_runtime.log"
open(runtime_log, "w") do io
    println(io, "timestamp,t,dt,umax")
end

# Zellzentrierte Geschwindigkeits-Magnitude
mag(I, u) = sqrt(sum(ntuple(i -> 0.25 * (u[I, i] + u[I+δ(i, I), i])^2, length(I))))

duration = 20.0
step = 0.1

anim = Plots.@animate for t in step:step:duration
    WaterLily.sim_step!(sim, t; remeasure=false)

    @inside sim.flow.σ[I] = mag(I, sim.flow.u)
    field = Array(sim.flow.σ[inside(sim.flow.σ)])

    Plots.heatmap(field';
        aspect_ratio=:equal,
        color=:viridis,
        clims=(0, 1),
        title="|u|, t=$(round(t, digits=2))",
        xlabel="x",
        ylabel="y"
    )

    dt = sim.flow.Δt[end]
    umax = maximum(Array(sim.flow.u[:, :, 1]))
    line = "$(Dates.format(Dates.now(), "HH:MM:SS")),$(round(t, digits=3)),$(round(dt, digits=6)),$(round(umax, digits=6))"

    open(runtime_log, "a") do io
        println(io, line)
    end
    println(line) # Live-Log im Terminal während der Berechnung
end

Plots.gif(anim, "$(case_id).gif"; fps=15)

# Residual-Log als Plot
plot_logger("$(case_id).log")
