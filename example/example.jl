import BasicModelInterface as BMI
using Heat
using TOML

path = normpath(@__DIR__, "../example/heat.toml")

m = BMI.initialize(Heat.Model, path)

BMI.update(m)
BMI.get_current_time(m)

temperature = BMI.get_value_ptr(m, "plate_surface__temperature")
