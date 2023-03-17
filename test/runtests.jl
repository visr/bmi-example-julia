import BasicModelInterface as BMI
using Heat
using Test
using TOML

@testset "Heat" begin

    @testset "initialize from defaults" begin
        m = BMI.initialize(Heat.Model)
        z0 = BMI.get_value_ptr(m, "plate_surface__temperature")
        @test all(z0 .>= 0)
        @test all(z0 .< 1)
    end

    @testset "initialize from dict" begin
        config = Dict("shape" => [7, 5])
        m = BMI.initialize(Heat.Model, config)
        ndim = BMI.get_grid_rank(m, 0)
        shape = zeros(Int, ndim)
        BMI.get_grid_shape(m, 0, shape)
        @test shape == [7, 5]
    end

    @testset "initialize from file" begin
        path = normpath(@__DIR__, "../example/heat.toml")
        config = TOML.parsefile(path)

        m = BMI.initialize(Heat.Model, path)
        ndim = BMI.get_grid_rank(m, 0)
        shape = zeros(Int, ndim)
        BMI.get_grid_shape(m, 0, shape)
        @test shape == config["shape"]
    end

    @testset "update" begin
        m = Heat.Model()
        n = 10
        for _ = 1:n
            BMI.update(m)
        end
        @test BMI.get_current_time(m) ≈ n * BMI.get_time_step(m)
    end

    @testset "update_until" begin
        m = Heat.Model()
        BMI.update_until(m, 10.1)
        @test BMI.get_current_time(m) ≈ 10.1
    end

    @testset "value copy" begin
        m = Heat.Model()

        dest0 = zeros(BMI.get_grid_size(m, 0))
        dest1 = zeros(BMI.get_grid_size(m, 0))

        z0 = BMI.get_value(m, "plate_surface__temperature", dest0)
        z1 = BMI.get_value(m, "plate_surface__temperature", dest1)

        @test z0 !== z1
        @test z0 == z1
    end

end  # testset "Heat"
