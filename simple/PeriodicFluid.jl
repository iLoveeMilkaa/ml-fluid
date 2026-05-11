import WaterLily
import Plots
import FFTW
import Random

# ---------------------------------------------------------------------------
# Smooth random velocity field via inverse FFT
#
# A random stream function ψ is built in Fourier space keeping only
# wavenumbers with |k| ≤ k_max (low frequencies → spatially smooth field).
# Spectral differentiation then yields a divergence-free velocity:
#   u =  ∂ψ/∂y  →  û[kx,ky] =  i·ky·(2π/L)·ψ̂[kx,ky]
#   v = -∂ψ/∂x  →  v̂[kx,ky] = -i·kx·(2π/L)·ψ̂[kx,ky]
# irfft maps back to physical space.
# ---------------------------------------------------------------------------
function random_velocity_field(L::Int; k_max::Int=4, seed::Int=42)
    rng = Random.MersenneTwister(seed)

    # rfft layout for (L,L) real input → output shape (L÷2+1, L)
    # dim-1 index kxp → kx = kxp-1   (kx ≥ 0, rfft conjugate symmetry)
    # dim-2 index jy  → ky = jy-1 for jy ≤ L÷2+1, else jy-1-L
    ψ̂ = zeros(ComplexF64, L÷2+1, L)
    for jy in 1:L
        ky = jy <= L÷2+1 ? jy-1 : jy-1-L
        for kx in 0:L÷2
            k2 = kx^2 + ky^2
            (k2 == 0 || k2 > k_max^2) && continue
            # amplitude ∝ 1/|k| → red spectrum, spatially smooth
            ψ̂[kx+1, jy] = (randn(rng) + im*randn(rng)) / sqrt(Float64(k2))
        end
    end

    # Spectral differentiation → divergence-free (u, v)
    κ  = 2π / L
    ûx = zeros(ComplexF64, L÷2+1, L)
    ûy = zeros(ComplexF64, L÷2+1, L)
    for jy in 1:L
        ky = jy <= L÷2+1 ? jy-1 : jy-1-L
        for kxp in 1:L÷2+1
            kx = kxp - 1
            ûx[kxp, jy] =  im * ky * κ * ψ̂[kxp, jy]
            ûy[kxp, jy] = -im * kx * κ * ψ̂[kxp, jy]
        end
    end

    ux = FFTW.irfft(ûx, L)
    uy = FFTW.irfft(ûy, L)

    U_max = maximum(sqrt.(ux .^ 2 .+ uy .^ 2))
    return ux ./ U_max, uy ./ U_max
end

function make_sim(; L=64, Re=1000.0, T=Float64, k_max=4)
    ν = T(L / Re)   # ν = U·L/Re with U=1

    ux, uy = random_velocity_field(L; k_max)

    # uBC initialises the field; periodic BCs take over afterwards
    uBC = (i, x, t) -> begin
        ix = mod1(round(Int, x[1]), L)
        iy = mod1(round(Int, x[2]), L)
        T(i == 1 ? ux[ix, iy] : uy[ix, iy])
    end

    WaterLily.Simulation((L, L), uBC, L; U=1, ν, T, perdir=(1, 2))
end

# Build simulation
sim = make_sim(L=64, Re=1000.0)

function umag(sim)  # cell-centre velocity magnitude
    u = sim.flow.u
    n = size(u, 1) - 2
    [
        sqrt(
            (0.5 * (u[i, j, 1] + u[i+1, j, 1]))^2 +
            (0.5 * (u[i, j, 2] + u[i, j+1, 2]))^2
        )
        for i in 2:n+1, j in 2:n+1
    ]
end

# Animated plot evolving in time
t_end = 20.0
step  = 0.5

anim = Plots.@animate for t in range(0.0, t_end; step)
    WaterLily.sim_step!(sim, t; remeasure=false)
    Plots.heatmap(
        umag(sim)';
        color=:viridis,
        aspect_ratio=:equal,
        clims=(0, 1),
        title="Absolute velocity |u|  t=$(round(WaterLily.sim_time(sim), digits=2))",
        xlabel="x",
        ylabel="y",
        colorbar_title="|u|",
    )
end

Plots.gif(anim, "plots/periodic_fluid.gif"; fps=10)
