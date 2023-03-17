"The 2D heat model."
module Heat

using Statistics
using TOML
using DSP: conv

"""
    struct Model

Solve the Heat equation on a grid.

# Keywords
- `shape::Tuple{Int, Int} = (10, 20)`
    The shape of the solution grid as (*columns*, *rows*).
- `spacing::Tuple{Float64, Float64} = (1.0, 1.0)`
    Spacing of grid columns and rows.
- `origin::Tuple{Float64, Float64} = (0.0, 0.0)`
    Coordinates of lower left corner of grid.
- `alpha::Float64 = 1.0`
    Thermal diffusivity.
- `time::Float64 = 0.0`
    Current model time.
- `time_step::Float64 = 0.0`
    Time step.
- `temperature::Matrix{Float64} = rand(shape...)`
    Temperature on the plate.
"""
Base.@kwdef mutable struct Model
    shape::Tuple{Int,Int} = (10, 20)
    spacing::Tuple{Float64,Float64} = (1.0, 1.0)
    origin::Tuple{Float64,Float64} = (0.0, 0.0)
    alpha::Float64 = 1.0
    time::Float64 = 0.0
    time_step::Float64 = minimum(spacing)^2 / (4.0 * alpha)
    temperature::Matrix{Float64} = rand(shape...)
end

function Base.show(io::IO, model::Model)
    time = model.time
    temp = mean(model.temperature)
    println(io, "Heat.Model at time = ", time, " and mean temperature = ", temp)
end

# convert the vectors from the TOML to tuples for the Model
to_tuple(x::AbstractVector) = Tuple(x)
to_tuple(x) = x

function Model(config::AbstractDict)
    kwargs = NamedTuple([Symbol(k) => to_tuple(v) for (k, v) in config])
    return Model(; kwargs...)
end

function Model(path::AbstractString)
    config = TOML.parsefile(path)
    return Model(config)
end

"Solve the 2D Heat Equation on a uniform mesh."
function solve_2d!(model::Model, time_step)
    (; temperature, spacing, alpha) = model

    dx2, dy2 = spacing .^ 2
    stencil =
        [
            0.0 dy2 0.0
            dx2 -2(dx2+dy2) dx2
            0.0 dy2 0.0
        ] .* alpha .* time_step ./ (2(dx2 * dy2))

    out = conv(temperature, stencil)
    cutout = @view out[begin+1:end-1, begin+1:end-1]
    temperature .+= cutout
    return model
end

function advance_in_time!(model::Model)
    solve_2d!(model, model.time_step)
    model.time += model.time_step
    return model
end

function advance_in_time!(model::Model, time_step)
    solve_2d!(model, time_step)
    model.time += time_step
    return model
end

include("bmi.jl")

end # module Heat
