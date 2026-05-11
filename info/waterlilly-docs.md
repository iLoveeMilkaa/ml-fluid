
<a id='WaterLily'></a>

<a id='WaterLily-1'></a>

# WaterLily


<a id='Introduction-and-Quickstart'></a>

<a id='Introduction-and-Quickstart-1'></a>

## Introduction and Quickstart

<a id='WaterLily' href='#WaterLily'>#</a>
**`WaterLily`** &mdash; *Module*.



**WaterLily.jl**

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://WaterLily-jl.github.io/WaterLily.jl/dev/) [![Examples](https://img.shields.io/badge/view-examples-blue.svg)](https://github.com/WaterLily-jl/WaterLily-Examples/) [![CI](https://github.com/WaterLily-jl/WaterLily.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/WaterLily-jl/WaterLily.jl/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/WaterLily-jl/WaterLily.jl/branch/master/graph/badge.svg?token=8XYFWKOUFN)](https://codecov.io/gh/WaterLily-jl/WaterLily.jl) [![DOI](https://zenodo.org/badge/DOI/10.1016/j.cpc.2025.109748.svg)](https://doi.org/10.1016/j.cpc.2025.109748)

![Julia flow](assets/julia.gif)

**Overview**

**WaterLily.jl** is a simple and fast fluid simulator written in pure Julia. This project is supported by awesome libraries developed within the Julia scientific community, and it aims to accelerate and enhance fluid simulations. Watch the JuliaCon2024 talk here:

[![JuliaCon2024 still and link](assets/JuliaCon2024.png)](https://www.youtube.com/watch?v=FwMh2rq9kOU)

If you have used WaterLily for research, please **cite us**! The [2025 paper](https://doi.org/10.1016/j.cpc.2025.109748) describes the main features of the solver and provides benchmarking, validation, and profiling results.

```julia
@article{WeymouthFont2025,
    author = {G.D. Weymouth and B. Font},
    title = {WaterLily.jl: A differentiable and backend-agnostic Julia solver for incompressible viscous flow around dynamic bodies},
    doi = {10.1016/j.cpc.2025.109748},
    journal = {Computer Physics Communications},
    year = {2025},
    volume = {315},
    pages = {109748},
}
```

**Method/capabilities**

WaterLily solves the unsteady incompressible 2D or 3D [Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations) on a Cartesian grid. The pressure Poisson equation is solved with a [geometric multigrid](https://en.wikipedia.org/wiki/Multigrid_method) method. Solid boundaries are modelled using the [Boundary Data Immersion Method](https://eprints.soton.ac.uk/369635/). The solver can run on serial CPU, multi-threaded CPU, or GPU backends.

**Example: Flow over a circle**

WaterLily lets the user can set the domain size and boundary conditions, the fluid viscosity (which determines the [Reynolds number](https://en.wikipedia.org/wiki/Reynolds_number)), and immerse solid obstacles. A large selection of examples, notebooks, and tutorials are found in the [WaterLily-Examples](https://github.com/WaterLily-jl/WaterLily-Examples) repository. Here, we will illustrate the basics by simulating and plotting the flow over a circle.

We define the size of the simulation domain as `n` by `m` cells. The circle has radius `m/8` and is centered at `(m/2,m/2)`. The flow boundary conditions are `(U,0)`, where we set `U=1`, and the Reynolds number is `Re=U*radius/ν` where `ν` (Greek "nu" U+03BD, not Latin lowercase "v") is the kinematic viscosity of the fluid.

```julia
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
```

The circle geometry is defined using a [signed distance function](https://en.wikipedia.org/wiki/Signed_distance_function#Applications). The `AutoBody` function uses [automatic differentiation](https://github.com/JuliaDiff/) to infer the other geometric parameters of the body automatically. Replace the circle's distance function with any other, and now you have the flow around something else... such as a [donut](https://github.com/WaterLily-jl/WaterLily-Examples/blob/main/examples/ThreeD_Donut.jl) or the [Julia logo](https://github.com/WaterLily-jl/WaterLily-Examples/blob/main/examples/TwoD_Julia.jl). For more complex geometries, [ParametricBodies.jl](https://github.com/WaterLily-jl/ParametricBodies.jl) defines a `body` using any parametric curve, such as a spline. See that repo (and the video above) for examples.

The code block above return a `Simulation` with the parameters we've defined. Now we can initialize a simulation (first line) and step it forward in time (second line)

```julia
circ = circle(3*2^5,2^6)
sim_step!(circ)
```

Note we've set `n,m` to be multiples of powers of 2, which is important when using the (very fast) geometric multi-grid solver.

We can now access and plot whatever variables we like. For example, we can plot the x-component of the velocity field using

```julia
using Plots
u = circ.flow.u[:,:,1] # first component is x
contourf(u') # transpose the array for the plot
```

![Initial velocity field](assets/u0.png)

As you can see, the velocity within the circle is zero, the velocity far from the circle is one, and there are accelerated and decelerated regions around the circle. The `sim_step!` has only taken a single time step, and this initial flow around our circle looks similar to the potential flow because the viscous boundary layer has not separated yet.

A set of [flow metric functions](https://github.com/WaterLily-jl/WaterLily.jl/blob/master/src/Metrics.jl) have been implemented, and we can use them to measure the simulation. The following code block defines a function to step the simulation to time `t` and then use the `pressure_force` metric to measure the force on the immersed body. The function is applied over a time range, and the forces are plotted.

```Julia
function get_forces!(sim,t)
    sim_step!(sim,t,remeasure=false)
    force = WaterLily.pressure_force(sim)
    force./(0.5sim.L*sim.U^2) # scale the forces!
end

# Simulate through the time range and get forces
time = 1:0.1:50 # time scale is sim.L/sim.U
forces = [get_forces!(circ,t) for t in time];

#Plot it
plot(time,[first.(forces) last.(forces)],
    labels=["drag" "lift"],
    xlabel="tU/L",
    ylabel="Pressure force coefficients")
```

![Pressure forces](assets/forces.png)

We can also plot the vorticity field instead of the u-velocity to see a snap-shot of the wake.

```julia
# Use curl(velocity) to compute vorticity `inside` the domain
ω = zeros(size(u));
@inside ω[I] = WaterLily.curl(3,I,circ.flow.u)*circ.L/circ.U

# Plot it using WaterLily's Plots Extension
flood(ω,clims = (-10,10),border=:none)
```

![Vorticity field](assets/vort.png)

Note that `flood` is a convience function within WaterLily to create 2D flood plots. As you can see, WaterLily correctly predicts that the flow is unsteady, with an alternating vortex street wake, leading to an oscillating side force and drag force.

**Multi-threading and GPU backends**

WaterLily uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) to multi-thread on CPU and run on GPU backends. The implementation method and speed-up are documented in the [2024 paper](https://physics.paperswithcode.com/paper/waterlily-jl-a-differentiable-and-backend), with costs as low as 1.44 nano-seconds measured per degree of freedom and time step!

Note that multi-threading requires *starting* Julia with the `--threads` argument, see [the multi-threading section](https://docs.julialang.org/en/v1/manual/multi-threading/) of the manual. If you are running Julia with multiple threads, KernelAbstractions will detect this and multi-thread the loops automatically.

Running on a GPU requires initializing the `Simulation` memory on the GPU, and care needs to be taken to move the data back to the CPU for visualization. As an example, let's compare a **3D** GPU simulation of a sphere to the **2D** multi-threaded CPU circle defined above

```Julia
using CUDA,WaterLily
function sphere(n,m;Re=100,U=1,T=Float64,mem=Array)
    radius, center = m/8, m/2-1
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)
    Simulation((n,m,m),(U,0,0), # 3D array size and BCs
                2radius;ν=U*2radius/Re,body, # no change
                T,   # Floating point type
                mem) # memory type
end

@assert CUDA.functional()      # is your CUDA GPU working??
GPUsim = sphere(3*2^5,2^6;T=Float32,mem=CuArray); # 3D GPU sim!
println(length(GPUsim.flow.u)) # 1.3M degrees-of freedom!
sim_step!(GPUsim)              # compile GPU code & run one step
@time sim_step!(GPUsim,50,remeasure=false) # 40s!!

CPUsim = circle(3*2^5,2^6);    # 2D CPU sim
println(length(CPUsim.flow.u)) # 0.013M degrees-of freedom!
sim_step!(CPUsim)              # compile GPU code & run one step
println(Threads.nthreads())    # I'm using 8 threads
@time sim_step!(CPUsim,50,remeasure=false) # 28s!!
```

As you can see, the 3D sphere set-up is almost identical to the 2D circle, but using 3D arrays means there are almost 1.3M degrees-of-freedom, 100x bigger than in 2D. Never the less, the simulation is quite fast on the GPU, only around 40% slower than the much smaller 2D simulation on a CPU with 8 threads. See the [2024 paper](https://physics.paperswithcode.com/paper/waterlily-jl-a-differentiable-and-backend) and the [examples repo](https://github.com/WaterLily-jl/WaterLily-Examples) for many more non-trivial examples including running on AMD GPUs.

Finally, KernelAbstractions does incur some CPU allocations for every loop, but other than this `sim_step!` is completely non-allocating. This is one reason why the speed-up improves as the size of the simulation increases.

**Contributing and issues**

We always appreciate new contributions, so please [submit a pull request](https://github.com/WaterLily-jl/WaterLily.jl/compare) with your changes and help us make WaterLily even better! Note that contributions need to be submitted together with benchmark results - WaterLily should always be fast! 😃 For this, we have a [fully automated benchmarking suite](https://github.com/WaterLily-jl/WaterLily-Benchmarks) that conducts performance tests. In short, to compare your changes with the latest WaterLily, clone the that repo and run the benchmarks with

```sh
git clone https://github.com/WaterLily-jl/WaterLily-Benchmarks && cd WaterLily-Benchmarks
sh benchmark.sh -wd "<your/waterlily/path>" -w "<your_waterlily_branch> master"
julia --project compare.jl
```

This will run benchmarks for CPU and GPU backends. If you do not have a GPU, simply pass `-b "Array"` when runnning `benchmark.sh`. More information on the benchmark suite is available in that [README](https://github.com/WaterLily-jl/WaterLily-Benchmarks/blob/main/README.md).

Of course, ideas, suggestions, and questions are welcome too! Please [raise an issue](https://github.com/WaterLily-jl/WaterLily.jl/issues/new/choose) to address any of these.

**Development goals**

  * Immerse obstacles defined by 3D meshes ([Meshing.jl](https://github.com/JuliaGeometry/Meshing.jl))
  * Multi-CPU/GPU simulations (https://github.com/WaterLily-jl/WaterLily.jl/pull/141)
  * Free-surface physics with ([Volume-of-Fluid](https://github.com/TzuYaoHuang/WaterLily.jl/blob/master/src/Multiphase.jl)) or other methods.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L1' class='documenter-source'>source</a><br>


<a id='Types-Methods-and-Functions'></a>

<a id='Types-Methods-and-Functions-1'></a>

## Types Methods and Functions



- [`WaterLily`](index.md#WaterLily)
- [`WaterLily.AbstractBody`](index.md#WaterLily.AbstractBody)
- [`WaterLily.AutoBody`](index.md#WaterLily.AutoBody)
- [`WaterLily.Flow`](index.md#WaterLily.Flow)
- [`WaterLily.MeanFlow`](index.md#WaterLily.MeanFlow)
- [`WaterLily.MultiLevelPoisson`](index.md#WaterLily.MultiLevelPoisson)
- [`WaterLily.NoBody`](index.md#WaterLily.NoBody)
- [`WaterLily.Poisson`](index.md#WaterLily.Poisson)
- [`WaterLily.RigidMap`](index.md#WaterLily.RigidMap)
- [`WaterLily.SetBody`](index.md#WaterLily.SetBody)
- [`WaterLily.Simulation`](index.md#WaterLily.Simulation)
- [`WaterLily.BC!`](index.md#WaterLily.BC!)
- [`WaterLily.CIj`](index.md#WaterLily.CIj-Union{Tuple{d}, Tuple{Any, CartesianIndex{d}, Any}} where d)
- [`WaterLily.GaussSeidelRB!`](index.md#WaterLily.GaussSeidelRB!-Union{Tuple{Poisson{T, S, V} where {S<:(AbstractArray{T}), V<:(AbstractArray{T})}}, Tuple{T}} where T)
- [`WaterLily.Jacobi!`](index.md#WaterLily.Jacobi!-Tuple{Any})
- [`WaterLily.L₂`](index.md#WaterLily.L₂-Tuple{Any})
- [`WaterLily.S`](index.md#WaterLily.S-Tuple{CartesianIndex{2}, Any})
- [`WaterLily._interp_clamp`](index.md#WaterLily._interp_clamp-Union{Tuple{T}, Tuple{D}, Tuple{StaticArraysCore.SVector{D, T}, NTuple{D, Int64}}} where {D, T})
- [`WaterLily.accelerate!`](index.md#WaterLily.accelerate!-Tuple{Any, Any, Nothing, Union{Nothing, Tuple}})
- [`WaterLily.apply!`](index.md#WaterLily.apply!-Tuple{Any, Any})
- [`WaterLily.check_nthreads`](index.md#WaterLily.check_nthreads-Tuple{})
- [`WaterLily.curl`](index.md#WaterLily.curl-Tuple{Any, Any, Any})
- [`WaterLily.curvature`](index.md#WaterLily.curvature-Tuple{AbstractMatrix})
- [`WaterLily.exitBC!`](index.md#WaterLily.exitBC!-Tuple{Any, Any, Any})
- [`WaterLily.inside`](index.md#WaterLily.inside-Tuple{AbstractArray})
- [`WaterLily.inside_u`](index.md#WaterLily.inside_u-Union{Tuple{N}, Tuple{NTuple{N, T} where T, Any}} where N)
- [`WaterLily.ke`](index.md#WaterLily.ke-Union{Tuple{m}, Tuple{CartesianIndex{m}, Any}, Tuple{CartesianIndex{m}, Any, Any}} where m)
- [`WaterLily.loc`](index.md#WaterLily.loc-Union{Tuple{N}, Tuple{Any, CartesianIndex{N}}, Tuple{Any, CartesianIndex{N}, Any}} where N)
- [`WaterLily.logger`](index.md#WaterLily.logger)
- [`WaterLily.measure`](index.md#WaterLily.measure-Tuple{AutoBody, Any, Any})
- [`WaterLily.measure!`](index.md#WaterLily.measure!-Union{Tuple{T}, Tuple{N}, Tuple{AbstractFlow{N, T}, AbstractBody}} where {N, T})
- [`WaterLily.measure!`](index.md#WaterLily.measure!)
- [`WaterLily.measure_sdf!`](index.md#WaterLily.measure_sdf!-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractBody}, Tuple{AbstractArray{T}, AbstractBody, Any}} where T)
- [`WaterLily.mom_correct!`](index.md#WaterLily.mom_correct!-Tuple{AbstractFlow, Any})
- [`WaterLily.mom_predict!`](index.md#WaterLily.mom_predict!-Tuple{AbstractFlow, Any, Any})
- [`WaterLily.mom_project!`](index.md#WaterLily.mom_project!-Tuple{AbstractFlow, AbstractPoisson, Any, Any})
- [`WaterLily.mom_step!`](index.md#WaterLily.mom_step!-Tuple{AbstractFlow, AbstractPoisson})
- [`WaterLily.mult!`](index.md#WaterLily.mult!-Tuple{Poisson, Any})
- [`WaterLily.nds`](index.md#WaterLily.nds-Union{Tuple{T}, Tuple{Any, AbstractVector{T}, Any}} where T)
- [`WaterLily.pcg!`](index.md#WaterLily.pcg!-Union{Tuple{Poisson{T, S, V} where {S<:(AbstractArray{T}), V<:(AbstractArray{T})}}, Tuple{T}} where T)
- [`WaterLily.perBC!`](index.md#WaterLily.perBC!-Tuple{Any, Tuple{}})
- [`WaterLily.perdot`](index.md#WaterLily.perdot-Tuple{Any, Any, Tuple{}})
- [`WaterLily.perturb!`](index.md#WaterLily.perturb!-Tuple{AbstractSimulation})
- [`WaterLily.pressure_force`](index.md#WaterLily.pressure_force-Tuple{Any})
- [`WaterLily.pressure_moment`](index.md#WaterLily.pressure_moment-Tuple{Any, Any})
- [`WaterLily.residual!`](index.md#WaterLily.residual!-Tuple{Poisson})
- [`WaterLily.sdf`](index.md#WaterLily.sdf)
- [`WaterLily.sdf`](index.md#WaterLily.sdf)
- [`WaterLily.sgs!`](index.md#WaterLily.sgs!-Tuple{Any, Any})
- [`WaterLily.sim_info`](index.md#WaterLily.sim_info-Tuple{AbstractSimulation})
- [`WaterLily.sim_step!`](index.md#WaterLily.sim_step!-Tuple{AbstractSimulation, Any})
- [`WaterLily.sim_time`](index.md#WaterLily.sim_time-Tuple{AbstractSimulation})
- [`WaterLily.slice`](index.md#WaterLily.slice-Union{Tuple{N}, Tuple{NTuple{N, T} where T, Any, Any}, Tuple{NTuple{N, T} where T, Any, Any, Any}} where N)
- [`WaterLily.solver!`](index.md#WaterLily.solver!-Tuple{Poisson})
- [`WaterLily.time`](index.md#WaterLily.time-Tuple{AbstractFlow})
- [`WaterLily.total_force`](index.md#WaterLily.total_force-Tuple{Any})
- [`WaterLily.udf!`](index.md#WaterLily.udf!-Tuple{Any, Nothing, Any})
- [`WaterLily.viscous_force`](index.md#WaterLily.viscous_force-Tuple{Any})
- [`WaterLily.δ`](index.md#WaterLily.δ-Union{Tuple{N}, Tuple{Any, Val{N}}} where N)
- [`WaterLily.λ₂`](index.md#WaterLily.λ₂-Tuple{CartesianIndex{3}, Any})
- [`WaterLily.ω`](index.md#WaterLily.ω-Tuple{CartesianIndex{3}, Any})
- [`WaterLily.ω_mag`](index.md#WaterLily.ω_mag-Tuple{CartesianIndex{3}, Any})
- [`WaterLily.ω_θ`](index.md#WaterLily.ω_θ-Tuple{CartesianIndex{3}, Any, Any, Any})
- [`WaterLily.∂`](index.md#WaterLily.∂-NTuple{4, Any})
- [`WaterLily.@inside`](index.md#WaterLily.@inside-Tuple{Any})
- [`WaterLily.@loop`](index.md#WaterLily.@loop-Tuple)

<a id='WaterLily.AbstractBody' href='#WaterLily.AbstractBody'>#</a>
**`WaterLily.AbstractBody`** &mdash; *Type*.



```julia
AbstractBody
```

Immersed body Abstract Type. Any `AbstractBody` subtype must implement

```julia
d,n,V = measure(body::AbstractBody, x, t=0, fastd²=Inf)
```

where `d` is the signed distance from `x` to the body at time `t`, and `n` & `V` are the normal and velocity vectors implied at `x`. A fast-approximate method can return `≈d,zero(x),zero(x)` if `d^2>fastd²`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Body.jl#L2-L12' class='documenter-source'>source</a><br>

<a id='WaterLily.AutoBody' href='#WaterLily.AutoBody'>#</a>
**`WaterLily.AutoBody`** &mdash; *Type*.



```julia
AutoBody(sdf,map=(x,t)->x) <: AbstractBody
```

  * `sdf(x::AbstractVector,t::Real)::Real`: signed distance function
  * `map(x::AbstractVector,t::Real)::AbstractVector`: coordinate mapping function

Implicitly define a geometry by its `sdf` and optional coordinate `map`. Note: the `map` is composed automatically i.e. `sdf(body::AutoBody,x,t) = body.sdf(body.map(x,t),t)`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/AutoBody.jl#L1-L9' class='documenter-source'>source</a><br>

<a id='WaterLily.Flow' href='#WaterLily.Flow'>#</a>
**`WaterLily.Flow`** &mdash; *Type*.



```julia
Flow{D::Int, T::Float, Sf<:AbstractArray{T,D}, Vf<:AbstractArray{T,D+1}, Tf<:AbstractArray{T,D+2}}
```

Composite type for a multidimensional immersed boundary flow simulation.

Flow solves the unsteady incompressible [Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations) on a Cartesian grid. Solid boundaries are modelled using the [Boundary Data Immersion Method](https://eprints.soton.ac.uk/369635/). The primary variables are the scalar pressure `p` (an array of dimension `D`) and the velocity vector field `u` (an array of dimension `D+1`).


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L76-L85' class='documenter-source'>source</a><br>

<a id='WaterLily.MeanFlow' href='#WaterLily.MeanFlow'>#</a>
**`WaterLily.MeanFlow`** &mdash; *Type*.



```julia
MeanFlow{T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Mf<:AbstractArray{T}}
```

Holds temporal averages of pressure, velocity, and squared-velocity tensor.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L154-L158' class='documenter-source'>source</a><br>

<a id='WaterLily.MultiLevelPoisson' href='#WaterLily.MultiLevelPoisson'>#</a>
**`WaterLily.MultiLevelPoisson`** &mdash; *Type*.



```julia
MultiLevelPoisson{N,M}
```

Composite type used to solve the pressure Poisson equation with a [geometric multigrid](https://en.wikipedia.org/wiki/Multigrid_method) method. The only variable is `levels`, a vector of nested `Poisson` systems.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/MultiLevelPoisson.jl#L38-L43' class='documenter-source'>source</a><br>

<a id='WaterLily.NoBody' href='#WaterLily.NoBody'>#</a>
**`WaterLily.NoBody`** &mdash; *Type*.



```julia
NoBody
```

Use for a simulation without a body.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Body.jl#L76-L80' class='documenter-source'>source</a><br>

<a id='WaterLily.Poisson' href='#WaterLily.Poisson'>#</a>
**`WaterLily.Poisson`** &mdash; *Type*.



```julia
Poisson{N,M}
```

Composite type for conservative variable coefficient Poisson equations:

```julia
∮ds β ∂x/∂n = σ
```

The resulting linear system is

```julia
Ax = [L+D+L']x = z
```

where A is symmetric, block-tridiagonal and extremely sparse. Moreover, `D[I]=-∑ᵢ(L[I,i]+L'[I,i])`. This means matrix storage, multiplication, ect can be easily implemented and optimized without external libraries.

To help iteratively solve the system above, the Poisson structure holds helper arrays for `inv(D)`, the error `ϵ`, and residual `r=z-Ax`. An iterative solution method then estimates the error `ϵ=̃A⁻¹r` and increments `x+=ϵ`, `r-=Aϵ`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L3-L21' class='documenter-source'>source</a><br>

<a id='WaterLily.RigidMap' href='#WaterLily.RigidMap'>#</a>
**`WaterLily.RigidMap`** &mdash; *Type*.



```julia
RigidMap(center, θ) <: AbstractBody
```

  * `x₀::SVector{D}`: coordinate of the center of the body
  * `θ::Union{Real, SVector{3}}`: rotation (single angle in 2D, and in 3D these are the rotation angle around                               the x, y, and z axes respectively.)
  * `V::SVector{D}=zero(center)`: linear velocity of the center
  * `xₚ::SVector{D}=zero(center)`: offset of the pivot point compared to center
  * `ω::Union{Real, SVector{3}}=zero(θ)`: angular velocity (scalar in 2D, vector in 3D)

Define a `RigidMap` for any `AbstractBody` using rigid body motion parameters.

RigidMap updates are computed externally via a set of ODEs and then updated in the simulation loop following:

```julia
using WaterLily,StaticArrays
body = AutoBody((x,t)->sqrt(sum(abs2,x))-4,RigidMap(SA{Float32}[16,16],0.f0;ω=0.1f0))
sim = Simulation((32,32),(1,0),8;body)
for n in 1:10
    # update body motion (example: constant angular velocity)
    θ = sim.body.map.θ + sim.body.map.ω*sim.flow.Δt[end]
    sim.body = setmap(sim.body; θ)
    # remeasure and step
    sim_step!(sim;remeasure=true)
end
```


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/RigidMap.jl#L1-L27' class='documenter-source'>source</a><br>

<a id='WaterLily.SetBody' href='#WaterLily.SetBody'>#</a>
**`WaterLily.SetBody`** &mdash; *Type*.



```julia
SetBody
```

Body defined as a lazy set operation on two `AbstractBody`s. The operations are only evaluated when `measure`d.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Body.jl#L85-L90' class='documenter-source'>source</a><br>

<a id='WaterLily.Simulation' href='#WaterLily.Simulation'>#</a>
**`WaterLily.Simulation`** &mdash; *Type*.



```julia
Simulation(dims::NTuple, uBC::Union{NTuple,Function}, L::Number;
           U=norm2(Uλ), Δt=0.25, ν=0., ϵ=1, g=nothing,
           perdir=(), exitBC=false,
           body::AbstractBody=NoBody(),
           T=Float32, mem=Array,
           flow_ctor=(dims,uBC;kw...)->Flow(dims,uBC;kw...),
           pois_ctor=flow->MultiLevelPoisson(flow.p,flow.μ₀,flow.σ;perdir))
```

Constructor for a WaterLily.jl simulation:

  * `dims`: Simulation domain dimensions.
  * `uBC`: Velocity field applied to boundary and acceleration conditions.     Define a `Tuple` for constant BCs, or a `Function` for space and time varying BCs `uBC(i,x,t)`.
  * `L`: Simulation length scale.
  * `U`: Simulation velocity scale. Required if using `Uλ::Function`.
  * `Δt`: Initial time step.
  * `ν`: Scaled viscosity (`Re=UL/ν`).
  * `g`: Domain acceleration, `g(i,x,t)=duᵢ/dt`
  * `ϵ`: BDIM kernel width.
  * `perdir`: Domain periodic boundary condition in the `(i,)` direction.
  * `uλ`: Velocity field applied to the initial condition.     Define a Tuple for homogeneous (per direction) IC, or a `Function` for space varying IC `uλ(i,x)`.
  * `exitBC`: Convective exit boundary condition in the `i=1` direction.
  * `body`: Immersed geometry.
  * `T`: Array element type.
  * `mem`: memory location. `Array`, `CuArray`, `ROCm` to run on CPU, NVIDIA, or AMD devices, respectively.
  * `flow_ctor`: Factory callable `(dims, uBC; kw...) -> AbstractFlow` to substitute a custom flow type.     The callable receives all standard keyword arguments forwarded from this constructor.     Used by downstream packages (e.g. LilyPad.jl) to inject a custom `AbstractFlow` subtype.
  * `pois_ctor`: Factory callable `flow -> AbstractPoisson` to substitute a custom Poisson solver.     Called after `flow_ctor` with the constructed flow as argument.     Used by downstream packages (e.g. BiotSavartBCs.jl) to inject a custom `AbstractPoisson` subtype.

See files in `examples` folder for examples.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L37-L72' class='documenter-source'>source</a><br>

<a id='WaterLily.BC!' href='#WaterLily.BC!'>#</a>
**`WaterLily.BC!`** &mdash; *Function*.



```julia
BC!(a,A)
```

Apply boundary conditions to the ghost cells of a *vector* field. A Dirichlet condition `a[I,i]=A[i]` is applied to the vector component *normal* to the domain boundary. For example `aₓ(x)=Aₓ ∀ x ∈ minmax(X)`. A zero Neumann condition is applied to the tangential components.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L208-L215' class='documenter-source'>source</a><br>

<a id='WaterLily.CIj-Union{Tuple{d}, Tuple{Any, CartesianIndex{d}, Any}} where d' href='#WaterLily.CIj-Union{Tuple{d}, Tuple{Any, CartesianIndex{d}, Any}} where d'>#</a>
**`WaterLily.CIj`** &mdash; *Method*.



```julia
CIj(j,I,k)
```

Replace jᵗʰ component of CartesianIndex with k


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L27-L30' class='documenter-source'>source</a><br>

<a id='WaterLily.GaussSeidelRB!-Union{Tuple{Poisson{T, S, V} where {S<:(AbstractArray{T}), V<:(AbstractArray{T})}}, Tuple{T}} where T' href='#WaterLily.GaussSeidelRB!-Union{Tuple{Poisson{T, S, V} where {S<:(AbstractArray{T}), V<:(AbstractArray{T})}}, Tuple{T}} where T'>#</a>
**`WaterLily.GaussSeidelRB!`** &mdash; *Method*.



```julia
GaussSeidelRB!(p::Poisson;it=4, ω=1)
```

Red-black Gauss-Seidel smoother. Runs `it` iterations; a complete red-black cycle requires `it` to be even. `ω` under-/over-relaxs the solution through scaling the deferred corrections in `increment!`. Note: This performs best on GPU configurations and is the default smoother.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L134-L140' class='documenter-source'>source</a><br>

<a id='WaterLily.Jacobi!-Tuple{Any}' href='#WaterLily.Jacobi!-Tuple{Any}'>#</a>
**`WaterLily.Jacobi!`** &mdash; *Method*.



```julia
Jacobi!(p::Poisson; it=1)
```

Jacobi smoother. Runs `it` iterations with relaxation parameter `ω` scaling the deferred corrections in `increment!`. Note: This runs for general backends but converges *very* slowly.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L105-L110' class='documenter-source'>source</a><br>

<a id='WaterLily.L₂-Tuple{Any}' href='#WaterLily.L₂-Tuple{Any}'>#</a>
**`WaterLily.L₂`** &mdash; *Method*.



```julia
L₂(a)
```

L₂ norm of array `a` excluding ghosts.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L63-L67' class='documenter-source'>source</a><br>

<a id='WaterLily.S-Tuple{CartesianIndex{2}, Any}' href='#WaterLily.S-Tuple{CartesianIndex{2}, Any}'>#</a>
**`WaterLily.S`** &mdash; *Method*.



```julia
S(I::CartesianIndex,u)
```

Rate-of-strain tensor.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L110-L114' class='documenter-source'>source</a><br>

<a id='WaterLily._interp_clamp-Union{Tuple{T}, Tuple{D}, Tuple{StaticArraysCore.SVector{D, T}, NTuple{D, Int64}}} where {D, T}' href='#WaterLily._interp_clamp-Union{Tuple{T}, Tuple{D}, Tuple{StaticArraysCore.SVector{D, T}, NTuple{D, Int64}}} where {D, T}'>#</a>
**`WaterLily._interp_clamp`** &mdash; *Method*.



```julia
interp(x::SVector, arr::AbstractArray)
```

Linear interpolation from array `arr` at Cartesian-coordinate `x`. Interpolation queries are clamped to the computational domain. Note: This routine works for any number of dimensions.

To interpolate from an `arr<:GPUArray`, the call for `interp` should be broadcasted over the coordinates `x` as follows:

```julia
p = CUDA.rand(10,18)
u = CUDA.rand(10,18,2)
x = CuArray([SA_F32[i-1.5, 2i+0.5] for i in 1:8])
WaterLily.interp.(x, Ref(p)) # Broadcast
WaterLily.interp.(x, Ref(u)) # Broadcast (x=[-0.5,2.5] is shifted to [0,2.5] because we are in a vector field)
```


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L271-L285' class='documenter-source'>source</a><br>

<a id='WaterLily.accelerate!-Tuple{Any, Any, Nothing, Union{Nothing, Tuple}}' href='#WaterLily.accelerate!-Tuple{Any, Any, Nothing, Union{Nothing, Tuple}}'>#</a>
**`WaterLily.accelerate!`** &mdash; *Method*.



```julia
accelerate!(r,t,g,U)
```

Accounts for applied and reference-frame acceleration using `rᵢ += g(i,x,t)+dU(i,x,t)/dt`


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L64-L68' class='documenter-source'>source</a><br>

<a id='WaterLily.apply!-Tuple{Any, Any}' href='#WaterLily.apply!-Tuple{Any, Any}'>#</a>
**`WaterLily.apply!`** &mdash; *Method*.



```julia
apply!(f, c)
```

Apply a vector function `f(i,x)` to the faces of a uniform staggered array `c` or a function `f(x)` to the center of a uniform array `c`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L188-L193' class='documenter-source'>source</a><br>

<a id='WaterLily.check_nthreads-Tuple{}' href='#WaterLily.check_nthreads-Tuple{}'>#</a>
**`WaterLily.check_nthreads`** &mdash; *Method*.



```julia
check_nthreads()
```

Check the number of threads available for the Julia session that loads WaterLily. A warning is shown when running in serial (JULIA*NUM*THREADS=1) with KernelAbstractions enabled.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L186-L191' class='documenter-source'>source</a><br>

<a id='WaterLily.curl-Tuple{Any, Any, Any}' href='#WaterLily.curl-Tuple{Any, Any, Any}'>#</a>
**`WaterLily.curl`** &mdash; *Method*.



```julia
curl(i,I,u)
```

Compute component `i` of $𝛁×𝐮$ at the **edge** of cell `I`. For example `curl(3,CartesianIndex(2,2,2),u)` will compute `ω₃(x=1.5,y=1.5,z=2)` as this edge produces the highest accuracy for this mix of cross derivatives on a staggered grid.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L53-L60' class='documenter-source'>source</a><br>

<a id='WaterLily.curvature-Tuple{AbstractMatrix}' href='#WaterLily.curvature-Tuple{AbstractMatrix}'>#</a>
**`WaterLily.curvature`** &mdash; *Method*.



```julia
curvature(A::AbstractMatrix)
```

Return `H,K` the mean and Gaussian curvature from `A=hessian(sdf)`. `K=tr(minor(A))` in 3D and `K=0` in 2D.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/AutoBody.jl#L40-L45' class='documenter-source'>source</a><br>

<a id='WaterLily.exitBC!-Tuple{Any, Any, Any}' href='#WaterLily.exitBC!-Tuple{Any, Any, Any}'>#</a>
**`WaterLily.exitBC!`** &mdash; *Method*.



```julia
exitBC!(u,u⁰,U,Δt)
```

Apply a 1D convection scheme to fill the ghost cell on the exit of the domain.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L237-L241' class='documenter-source'>source</a><br>

<a id='WaterLily.inside-Tuple{AbstractArray}' href='#WaterLily.inside-Tuple{AbstractArray}'>#</a>
**`WaterLily.inside`** &mdash; *Method*.



```julia
inside(a;buff=1)
```

Return CartesianIndices range excluding a single layer of cells on all boundaries.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L42-L46' class='documenter-source'>source</a><br>

<a id='WaterLily.inside_u-Union{Tuple{N}, Tuple{NTuple{N, T} where T, Any}} where N' href='#WaterLily.inside_u-Union{Tuple{N}, Tuple{NTuple{N, T} where T, Any}} where N'>#</a>
**`WaterLily.inside_u`** &mdash; *Method*.



```julia
inside_u(dims,j)
```

Return CartesianIndices range excluding the ghost-cells on the boundaries of a *vector* array on face `j` with size `dims`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L49-L54' class='documenter-source'>source</a><br>

<a id='WaterLily.ke-Union{Tuple{m}, Tuple{CartesianIndex{m}, Any}, Tuple{CartesianIndex{m}, Any, Any}} where m' href='#WaterLily.ke-Union{Tuple{m}, Tuple{CartesianIndex{m}, Any}, Tuple{CartesianIndex{m}, Any, Any}} where m'>#</a>
**`WaterLily.ke`** &mdash; *Method*.



```julia
ke(I::CartesianIndex,u,U=0)
```

Compute $½∥𝐮-𝐔∥²$ at center of cell `I` where `U` can be used to subtract a background flow (by default, `U=0`).


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L20-L25' class='documenter-source'>source</a><br>

<a id='WaterLily.loc-Union{Tuple{N}, Tuple{Any, CartesianIndex{N}}, Tuple{Any, CartesianIndex{N}, Any}} where N' href='#WaterLily.loc-Union{Tuple{N}, Tuple{Any, CartesianIndex{N}}, Tuple{Any, CartesianIndex{N}, Any}} where N'>#</a>
**`WaterLily.loc`** &mdash; *Method*.



```julia
loc(i,I) = loc(Ii)
```

Location in space of the cell at CartesianIndex `I` at face `i`. Using `i=0` returns the cell center s.t. `loc = I`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L178-L183' class='documenter-source'>source</a><br>

<a id='WaterLily.logger' href='#WaterLily.logger'>#</a>
**`WaterLily.logger`** &mdash; *Function*.



```julia
logger(fname="WaterLily")
```

Set up a logger to write the pressure solver data to a logging file named `WaterLily.log`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L11-L15' class='documenter-source'>source</a><br>

<a id='WaterLily.measure!' href='#WaterLily.measure!'>#</a>
**`WaterLily.measure!`** &mdash; *Function*.



```julia
measure!(sim::Simulation,t=timeNext(sim))
```

Measure a dynamic `body` to update the `flow` and `pois` coefficients.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L129-L133' class='documenter-source'>source</a><br>

<a id='WaterLily.measure!-Union{Tuple{T}, Tuple{N}, Tuple{AbstractFlow{N, T}, AbstractBody}} where {N, T}' href='#WaterLily.measure!-Union{Tuple{T}, Tuple{N}, Tuple{AbstractFlow{N, T}, AbstractBody}} where {N, T}'>#</a>
**`WaterLily.measure!`** &mdash; *Method*.



```julia
measure!(flow::AbstractFlow, body::AbstractBody; t=0, ϵ=1)
```

Queries the body geometry to fill the arrays:

  * `flow.μ₀`, Zeroth kernel moment
  * `flow.μ₁`, First kernel moment scaled by the body normal
  * `flow.V`,  Body velocity

at time `t` using an immersion kernel of size `ϵ`. The velocity is only filled within a narrow band of size `2+ϵ` around the body. This function also fills `flow.σ` with the signed distance function.

See Maertens & Weymouth, doi:[10.1016/j.cma.2014.09.007](https://doi.org/10.1016/j.cma.2014.09.007).


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Body.jl#L14-L27' class='documenter-source'>source</a><br>

<a id='WaterLily.measure-Tuple{AutoBody, Any, Any}' href='#WaterLily.measure-Tuple{AutoBody, Any, Any}'>#</a>
**`WaterLily.measure`** &mdash; *Method*.



```julia
d,n,V = measure(body::AutoBody,x,t;fastd²=Inf)
```

Determine the implicit geometric properties from the `sdf` and `map`. The gradient of `d=sdf(map(x,t))` is used to improve `d` for pseudo-sdfs. The velocity is determined *solely* from the optional `map` function. Skips the `n,V` calculation when `d²>fastd²`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/AutoBody.jl#L21-L28' class='documenter-source'>source</a><br>

<a id='WaterLily.measure_sdf!-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractBody}, Tuple{AbstractArray{T}, AbstractBody, Any}} where T' href='#WaterLily.measure_sdf!-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractBody}, Tuple{AbstractArray{T}, AbstractBody, Any}} where T'>#</a>
**`WaterLily.measure_sdf!`** &mdash; *Method*.



```julia
measure_sdf!(a::AbstractArray, body::AbstractBody, t=0; fastd²=0)
```

Uses `sdf(body,x,t)` to fill `a`. Defaults to fastd²=0 for quick evaluation.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Body.jl#L69-L73' class='documenter-source'>source</a><br>

<a id='WaterLily.mom_correct!-Tuple{AbstractFlow, Any}' href='#WaterLily.mom_correct!-Tuple{AbstractFlow, Any}'>#</a>
**`WaterLily.mom_correct!`** &mdash; *Method*.



```julia
mom_correct!(a::AbstractFlow, t; λ=quick, udf=nothing, kwargs...)
```

Corrector phase of `mom_step!`: advect under the projected `u`, apply BDIM, blend with the trapezoidal weight, enforce BCs at time-step end-time `t`. On return `a.u` is BC-consistent and ready for pressure projection.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L168-L174' class='documenter-source'>source</a><br>

<a id='WaterLily.mom_predict!-Tuple{AbstractFlow, Any, Any}' href='#WaterLily.mom_predict!-Tuple{AbstractFlow, Any, Any}'>#</a>
**`WaterLily.mom_predict!`** &mdash; *Method*.



```julia
mom_predict!(a::AbstractFlow, t₀, t₁; λ=quick, udf=nothing, kwargs...)
```

Predictor phase of `mom_step!`: advect under `u⁰`, apply BDIM, enforce BCs. On return `a.u` is BC-consistent and ready for pressure projection. `t₀` and `t₁` are the start and end times of the step; BCs are enforced at the end-of-step time `sum(a.Δt)`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L152-L159' class='documenter-source'>source</a><br>

<a id='WaterLily.mom_project!-Tuple{AbstractFlow, AbstractPoisson, Any, Any}' href='#WaterLily.mom_project!-Tuple{AbstractFlow, AbstractPoisson, Any, Any}'>#</a>
**`WaterLily.mom_project!`** &mdash; *Method*.



```julia
mom_project!(a::AbstractFlow, b::AbstractPoisson, w, t)
```

Projection phase of `mom_step!`: solve the pressure Poisson equation, correct the velocity by `w·Δt·∇p`, and re-enforce BCs. On return `a.u` is divergence-free and BC-consistent.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L183-L189' class='documenter-source'>source</a><br>

<a id='WaterLily.mom_step!-Tuple{AbstractFlow, AbstractPoisson}' href='#WaterLily.mom_step!-Tuple{AbstractFlow, AbstractPoisson}'>#</a>
**`WaterLily.mom_step!`** &mdash; *Method*.



```julia
mom_step!(a::AbstractFlow,b::AbstractPoisson;λ=quick,udf=nothing,kwargs...)
```

Integrate the `Flow` one time step using the [Boundary Data Immersion Method](https://eprints.soton.ac.uk/369635/) and the `AbstractPoisson` pressure solver to project the velocity onto an incompressible flow.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L120-L125' class='documenter-source'>source</a><br>

<a id='WaterLily.mult!-Tuple{Poisson, Any}' href='#WaterLily.mult!-Tuple{Poisson, Any}'>#</a>
**`WaterLily.mult!`** &mdash; *Method*.



```julia
mult!(p::Poisson,x)
```

Efficient function for Poisson matrix-vector multiplication. Fills `p.z = p.A x` with 0 in the ghost cells.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L57-L62' class='documenter-source'>source</a><br>

<a id='WaterLily.nds-Union{Tuple{T}, Tuple{Any, AbstractVector{T}, Any}} where T' href='#WaterLily.nds-Union{Tuple{T}, Tuple{Any, AbstractVector{T}, Any}} where T'>#</a>
**`WaterLily.nds`** &mdash; *Method*.



```julia
nds(body,x,t)
```

BDIM-masked surface normal.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L86-L90' class='documenter-source'>source</a><br>

<a id='WaterLily.pcg!-Union{Tuple{Poisson{T, S, V} where {S<:(AbstractArray{T}), V<:(AbstractArray{T})}}, Tuple{T}} where T' href='#WaterLily.pcg!-Union{Tuple{Poisson{T, S, V} where {S<:(AbstractArray{T}), V<:(AbstractArray{T})}}, Tuple{T}} where T'>#</a>
**`WaterLily.pcg!`** &mdash; *Method*.



```julia
pcg!(p::Poisson; it=6)
```

Conjugate-Gradient smoother with Jacobi predictioning. Runs at most `it` iterations, but will exit early if the Gram-Schmidt update parameter `|α| < 1%` or `|r D⁻¹ r| < 1e-8`. Note: This runs for general backends.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L151-L157' class='documenter-source'>source</a><br>

<a id='WaterLily.perBC!-Tuple{Any, Tuple{}}' href='#WaterLily.perBC!-Tuple{Any, Tuple{}}'>#</a>
**`WaterLily.perBC!`** &mdash; *Method*.



```julia
perBC!(a,perdir)
```

Apply periodic conditions to the ghost cells of a *scalar* field.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L250-L254' class='documenter-source'>source</a><br>

<a id='WaterLily.perdot-Tuple{Any, Any, Tuple{}}' href='#WaterLily.perdot-Tuple{Any, Any, Tuple{}}'>#</a>
**`WaterLily.perdot`** &mdash; *Method*.



```julia
perdot(a,b,perdir)
```

Apply dot product to the inner cells of two *scalar* fields, assuming zero values in ghost cell when using Neumann BC.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L262-L266' class='documenter-source'>source</a><br>

<a id='WaterLily.perturb!-Tuple{AbstractSimulation}' href='#WaterLily.perturb!-Tuple{AbstractSimulation}'>#</a>
**`WaterLily.perturb!`** &mdash; *Method*.



```julia
perturb!(sim; noise=0.1)
```

Perturb the velocity field of a simulation with `noise` level with respect to velocity scale `U`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L145-L148' class='documenter-source'>source</a><br>

<a id='WaterLily.pressure_force-Tuple{Any}' href='#WaterLily.pressure_force-Tuple{Any}'>#</a>
**`WaterLily.pressure_force`** &mdash; *Method*.



```julia
pressure_force(sim::Simulation)
```

Compute the pressure force on an immersed body.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L96-L100' class='documenter-source'>source</a><br>

<a id='WaterLily.pressure_moment-Tuple{Any, Any}' href='#WaterLily.pressure_moment-Tuple{Any, Any}'>#</a>
**`WaterLily.pressure_moment`** &mdash; *Method*.



```julia
pressure_moment(x₀,sim::Simulation)
```

Computes the pressure moment on an immersed body relative to point x₀.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L140-L144' class='documenter-source'>source</a><br>

<a id='WaterLily.residual!-Tuple{Poisson}' href='#WaterLily.residual!-Tuple{Poisson}'>#</a>
**`WaterLily.residual!`** &mdash; *Method*.



```julia
residual!(p::Poisson)
```

Computes the residual `r = z-Ax` and corrects it such that `r = 0` if `iD==0` which ensures local satisfiability     and `sum(r) = 0` which ensures global satisfiability.

The global correction is done by adjusting all points uniformly, minimizing the local effect. Other approaches are possible.

Note: These corrections mean `x` is not strictly solving `Ax=z`, but without the corrections, no solution exists.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L78-L91' class='documenter-source'>source</a><br>

<a id='WaterLily.sdf' href='#WaterLily.sdf'>#</a>
**`WaterLily.sdf`** &mdash; *Function*.



```julia
d = sdf(body::AutoBody,x,t) = body.sdf(body.map(x,t),t)
```


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/AutoBody.jl#L16-L18' class='documenter-source'>source</a><br>

<a id='WaterLily.sdf' href='#WaterLily.sdf'>#</a>
**`WaterLily.sdf`** &mdash; *Function*.



```julia
d = sdf(a::AbstractBody,x,t=0;fastd²=0)
```

Measure only the distance. Defaults to fastd²=0 for quick evaluation.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Body.jl#L62-L66' class='documenter-source'>source</a><br>

<a id='WaterLily.sgs!-Tuple{Any, Any}' href='#WaterLily.sgs!-Tuple{Any, Any}'>#</a>
**`WaterLily.sgs!`** &mdash; *Method*.



```julia
sgs!(flow, t; νₜ, S, Cs, Δ)
```

Implements a user-defined function `udf` to model subgrid-scale LES stresses based on the Boussinesq approximation     τᵃᵢⱼ = τʳᵢⱼ - (1/3)τʳₖₖδᵢⱼ = -2νₜS̅ᵢⱼ where             ▁▁▁▁     τʳᵢⱼ =  uᵢuⱼ - u̅ᵢu̅ⱼ

and we add -∂ⱼ(τᵃᵢⱼ) to the RHS as a body force (the isotropic part of the tensor is automatically modelled by the pressure gradient term). Users need to define the turbulent viscosity function `νₜ` and pass it as a keyword argument to this function together with rate-of-strain tensor array buffer `S`, Smagorinsky constant `Cs`, and filter width `Δ`. For example, the standard Smagorinsky–Lilly model for the sub-grid scale stresses is

```julia
νₜ = (CₛΔ)²|S̅ᵢⱼ|=(CₛΔ)²√(2S̅ᵢⱼS̅ᵢⱼ)
```

It can be implemented as     `smagorinsky(I::CartesianIndex{m} where m; S, Cs, Δ) = @views (Cs*Δ)^2*sqrt(dot(S[I,:,:],S[I,:,:]))` and passed into `sim_step!` as a keyword argument together with the varibles than the function needs (`S`, `Cs`, and `Δ`):     `sim_step!(sim, ...; udf=sgs, νₜ=smagorinsky, S, Cs, Δ)`


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L314-L334' class='documenter-source'>source</a><br>

<a id='WaterLily.sim_info-Tuple{AbstractSimulation}' href='#WaterLily.sim_info-Tuple{AbstractSimulation}'>#</a>
**`WaterLily.sim_info`** &mdash; *Method*.



```julia
sim_info(sim::AbstractSimulation)
```

Prints information on the current state of a simulation.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L139-L142' class='documenter-source'>source</a><br>

<a id='WaterLily.sim_step!-Tuple{AbstractSimulation, Any}' href='#WaterLily.sim_step!-Tuple{AbstractSimulation, Any}'>#</a>
**`WaterLily.sim_step!`** &mdash; *Method*.



```julia
sim_step!(sim::AbstractSimulation,t_end;remeasure=true,λ=quick,max_steps=typemax(Int),verbose=false,
    udf=nothing,kwargs...)
```

Integrate the simulation `sim` up to dimensionless time `t_end`. If `remeasure=true`, the body is remeasured at every time step. Can be set to `false` for static geometries to speed up simulation. A user-defined function `udf` can be passed to arbitrarily modify the `::AbstractFlow` during the predictor and corrector steps. If the `udf` user keyword arguments, these needs to be included in the `sim_step!` call as well. A `λ::Function` function can be passed as a custom convective scheme, following the interface of `λ(u,c,d)` (for upstream, central, downstream points).


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L105-L115' class='documenter-source'>source</a><br>

<a id='WaterLily.sim_time-Tuple{AbstractSimulation}' href='#WaterLily.sim_time-Tuple{AbstractSimulation}'>#</a>
**`WaterLily.sim_time`** &mdash; *Method*.



```julia
sim_time(sim::Simulation)
```

Return the current dimensionless time of the simulation `tU/L` where `t=sum(Δt)`, and `U`,`L` are the simulation velocity and length scales.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/WaterLily.jl#L96-L102' class='documenter-source'>source</a><br>

<a id='WaterLily.slice-Union{Tuple{N}, Tuple{NTuple{N, T} where T, Any, Any}, Tuple{NTuple{N, T} where T, Any, Any, Any}} where N' href='#WaterLily.slice-Union{Tuple{N}, Tuple{NTuple{N, T} where T, Any, Any}, Tuple{NTuple{N, T} where T, Any, Any, Any}} where N'>#</a>
**`WaterLily.slice`** &mdash; *Method*.



```julia
slice(dims,i,j,low=1)
```

Return `CartesianIndices` range slicing through an array of size `dims` in dimension `j` at index `i`. `low` optionally sets the lower extent of the range in the other dimensions.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L197-L203' class='documenter-source'>source</a><br>

<a id='WaterLily.solver!-Tuple{Poisson}' href='#WaterLily.solver!-Tuple{Poisson}'>#</a>
**`WaterLily.solver!`** &mdash; *Method*.



```julia
solver!(A::Poisson;tol=1e-4,itmx=1e3)
```

Approximate iterative solver for the Poisson matrix equation `Ax=b`.

  * `A`: Poisson matrix with working arrays.
  * `A.x`: Solution vector. Can start with an initial guess.
  * `A.z`: Right-Hand-Side vector. Will be overwritten!
  * `A.n[end]`: stores the number of iterations performed.
  * `tol`: Convergence tolerance on the `L₂`-norm residual.
  * `itmx`: Maximum number of iterations.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Poisson.jl#L183-L194' class='documenter-source'>source</a><br>

<a id='WaterLily.time-Tuple{AbstractFlow}' href='#WaterLily.time-Tuple{AbstractFlow}'>#</a>
**`WaterLily.time`** &mdash; *Method*.



```julia
time(a::AbstractFlow)
```

Current flow time.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L139-L143' class='documenter-source'>source</a><br>

<a id='WaterLily.total_force-Tuple{Any}' href='#WaterLily.total_force-Tuple{Any}'>#</a>
**`WaterLily.total_force`** &mdash; *Method*.



```julia
total_force(sim::Simulation)
```

Compute the total force on an immersed body.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L132-L136' class='documenter-source'>source</a><br>

<a id='WaterLily.udf!-Tuple{Any, Nothing, Any}' href='#WaterLily.udf!-Tuple{Any, Nothing, Any}'>#</a>
**`WaterLily.udf!`** &mdash; *Method*.



```julia
udf!(flow::AbstractFlow,udf::Function,t)
```

User defined function using `udf::Function` to operate on `flow::AbstractFlow` during the predictor and corrector step, in sync with time `t`. Keyword arguments must be passed to `sim_step!` for them to be carried over the actual function call.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Flow.jl#L213-L218' class='documenter-source'>source</a><br>

<a id='WaterLily.viscous_force-Tuple{Any}' href='#WaterLily.viscous_force-Tuple{Any}'>#</a>
**`WaterLily.viscous_force`** &mdash; *Method*.



```julia
viscous_force(sim::Simulation)
```

Compute the viscous force on an immersed body.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L118-L122' class='documenter-source'>source</a><br>

<a id='WaterLily.δ-Union{Tuple{N}, Tuple{Any, Val{N}}} where N' href='#WaterLily.δ-Union{Tuple{N}, Tuple{Any, Val{N}}} where N'>#</a>
**`WaterLily.δ`** &mdash; *Method*.



```julia
δ(i,N::Int)
δ(i,I::CartesianIndex{N}) where {N}
```

Return a CartesianIndex of dimension `N` which is one at index `i` and zero elsewhere.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L33-L38' class='documenter-source'>source</a><br>

<a id='WaterLily.λ₂-Tuple{CartesianIndex{3}, Any}' href='#WaterLily.λ₂-Tuple{CartesianIndex{3}, Any}'>#</a>
**`WaterLily.λ₂`** &mdash; *Method*.



```julia
λ₂(I::CartesianIndex{3},u)
```

λ₂ is a deformation tensor metric to identify vortex cores. See [https://en.wikipedia.org/wiki/Lambda2_method](https://en.wikipedia.org/wiki/Lambda2_method) and Jeong, J., & Hussain, F., doi:[10.1017/S0022112095000462](https://doi.org/10.1017/S0022112095000462)


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L40-L46' class='documenter-source'>source</a><br>

<a id='WaterLily.ω-Tuple{CartesianIndex{3}, Any}' href='#WaterLily.ω-Tuple{CartesianIndex{3}, Any}'>#</a>
**`WaterLily.ω`** &mdash; *Method*.



```julia
ω(I::CartesianIndex{3},u)
```

Compute 3-vector $𝛚=𝛁×𝐮$ at the center of cell `I`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L62-L66' class='documenter-source'>source</a><br>

<a id='WaterLily.ω_mag-Tuple{CartesianIndex{3}, Any}' href='#WaterLily.ω_mag-Tuple{CartesianIndex{3}, Any}'>#</a>
**`WaterLily.ω_mag`** &mdash; *Method*.



```julia
ω_mag(I::CartesianIndex{3},u)
```

Compute $∥𝛚∥$ at the center of cell `I`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L68-L72' class='documenter-source'>source</a><br>

<a id='WaterLily.ω_θ-Tuple{CartesianIndex{3}, Any, Any, Any}' href='#WaterLily.ω_θ-Tuple{CartesianIndex{3}, Any, Any, Any}'>#</a>
**`WaterLily.ω_θ`** &mdash; *Method*.



```julia
ω_θ(I::CartesianIndex{3},z,center,u)
```

Compute $𝛚⋅𝛉$ at the center of cell `I` where $𝛉$ is the azimuth direction around vector `z` passing through `center`.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L74-L79' class='documenter-source'>source</a><br>

<a id='WaterLily.∂-NTuple{4, Any}' href='#WaterLily.∂-NTuple{4, Any}'>#</a>
**`WaterLily.∂`** &mdash; *Method*.



```julia
∂(i,j,I,u)
```

Compute $∂uᵢ/∂xⱼ$ at center of cell `I`. Cross terms are computed less accurately than inline terms because of the staggered grid.


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/Metrics.jl#L29-L34' class='documenter-source'>source</a><br>

<a id='WaterLily.@inside-Tuple{Any}' href='#WaterLily.@inside-Tuple{Any}'>#</a>
**`WaterLily.@inside`** &mdash; *Macro*.



```julia
@inside <expr>
```

Simple macro to automate efficient loops over cells excluding ghosts. For example,

```julia
@inside p[I] = sum(loc(0,I))
```

becomes

```julia
@loop p[I] = sum(loc(0,I)) over I ∈ inside(p)
```

See [`@loop`](index.md#WaterLily.@loop-Tuple).


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L70-L82' class='documenter-source'>source</a><br>

<a id='WaterLily.@loop-Tuple' href='#WaterLily.@loop-Tuple'>#</a>
**`WaterLily.@loop`** &mdash; *Macro*.



```julia
@loop <expr> over <I ∈ R>
```

Macro to automate fast loops using @simd when running in serial, or KernelAbstractions when running multi-threaded CPU or GPU.

For example

```julia
@loop a[I,i] += sum(loc(i,I)) over I ∈ R
```

becomes

```julia
@simd for I ∈ R
    @fastmath @inbounds a[I,i] += sum(loc(i,I))
end
```

on serial execution, or

```julia
@kernel function kern(a,i,@Const(I0))
    I ∈ @index(Global,Cartesian)+I0
    @fastmath @inbounds a[I,i] += sum(loc(i,I))
end
kern(get_backend(a),64)(a,i,R[1]-oneunit(R[1]),ndrange=size(R))
```

when multi-threading on CPU or using CuArrays. Note that `get_backend` is used on the *first* variable in `expr` (`a` in this example).


<a target='_blank' href='https://github.com/WaterLily-jl/WaterLily.jl/blob/de13bc834fd8ffe166eb4eed3edc4dbd5e3bc0e9/src/util.jl#L105-L131' class='documenter-source'>source</a><br>

