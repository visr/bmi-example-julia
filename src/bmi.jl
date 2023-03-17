import BasicModelInterface as BMI

const varname = "plate_surface__temperature"

BMI.initialize(::Type{Model}) = Model()
BMI.initialize(::Type{Model}, x) = Model(x)

BMI.update(m::Model) = advance_in_time!(m)

function BMI.update_until(m::Model, t)
    time_step = BMI.get_time_step(m)
    n_steps = (t - BMI.get_current_time(m)) / time_step

    n_whole_steps = floor(Int, n_steps)
    fractional_step = rem(n_steps, 1) * time_step

    for _ = 1:n_whole_steps
        advance_in_time!(m)
    end

    if fractional_step > 0
        advance_in_time!(m, fractional_step)
    end
end

hasvar(name) = name == varname || throw(KeyError(name))
hasgrid(grid) = grid == 0 || throw(KeyError(grid))

BMI.finalize(m::Model) = m

BMI.get_var_type(m::Model, name) = hasvar(name) && repr(eltype(m.temperature))
BMI.get_var_units(m::Model, name) = hasvar(name) && "K"
BMI.get_var_nbytes(m::Model, name) = hasvar(name) && sizeof(m.temperature)
BMI.get_var_itemsize(m::Model, name) = hasvar(name) && sizeof(eltype(m.temperature))
BMI.get_var_location(m::Model, name) = hasvar(name) && "node"
BMI.get_var_grid(m::Model, name) = hasvar(name) && 0

BMI.get_grid_rank(m::Model, grid) = hasgrid(grid) && ndims(m.temperature)
BMI.get_grid_size(m::Model, grid) = hasgrid(grid) && length(m.temperature)

BMI.get_value_ptr(m::Model, name) = hasvar(name) && m.temperature

function BMI.get_value(m::Model, name, dest)
    val = BMI.get_value_ptr(m, name)
    copyto!(dest, val)
end

function BMI.get_value(m::Model, name, dest, inds)
    val = BMI.get_value_ptr(m, name)
    copyto!(dest, val[inds])
end

function BMI.set_value(m::Model, name, src)
    val = BMI.get_value_ptr(m, name)
    val .= src
    m
end

function BMI.set_value_at_indices(m::Model, name, inds, src)
    val = BMI.get_value_ptr(m, name)
    val[inds] .= src
    m
end

BMI.get_component_name(m::Model) = "The 2D Heat Equation"

BMI.get_input_item_count(m::Model) = 1
BMI.get_output_item_count(m::Model) = 1

BMI.get_input_var_names(m::Model) = [varname]
BMI.get_output_var_names(m::Model) = [varname]

function BMI.get_grid_shape(m::Model, grid, shape)
    hasgrid(grid)
    shape .= size(m.temperature)
end

function BMI.get_grid_spacing(m::Model, grid, spacing)
    hasgrid(grid)
    spacing .= m.spacing
end

function BMI.get_grid_origin(m::Model, grid, origin)
    hasgrid(grid)
    origin .= m.origin
end

BMI.get_grid_type(m::Model, grid) = hasgrid(grid) && "uniform_rectilinear"
BMI.get_start_time(m::Model) = 0.0
BMI.get_end_time(m::Model) = Inf
BMI.get_current_time(m::Model) = m.time
BMI.get_time_step(m::Model) = m.time_step
BMI.get_time_units(m::Model) = "s"
BMI.get_grid_node_count(m::Model, grid) = BMI.get_grid_size(m, grid)
