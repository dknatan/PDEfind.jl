using PDEfind
using Test

@testset "grid_operations.jl" begin
    time_array, space_array = collect(0.0:0.03:0.4), collect(0.0:0.1:4.3)
    nt, nx = length(time_array), length(space_array)
    g = TimeSpaceGrid1D(time_array, space_array)
    @test typeof(g) == TimeSpaceGrid1D{Float64}
    ∂t = TimeDerivative(g)
    @test size(∂t.M) == (nt, nt)
    ∂x = XDerivative(g)
    @test size(∂x.M) == (nx, nx)

    data_array = time_array * space_array'
    dt_data_array = ∂t(data_array) 
    dt_data_array_flat = reshape(dt_data_array, g.N) 
    max_derivative_degree = 4
    variables_matrix = GetVariablesMatrix(max_derivative_degree, data_array, g, ∂x)
    @test length(variables_matrix) == max_derivative_degree + 1
    @test length(variables_matrix[1]) == g.N
end    


@testset "function_bases.jl" begin
    max_derivative_degree = 2
    max_poly_degree = 3
    n_variables = max_derivative_degree + 1 # include 0th-derivative i.e. function itself
    MyBasis = PolynomialBasis(max_poly_degree, n_variables, skip_constant = false)
    @test typeof(MyBasis) == PolynomialBasis
    @test length(MyBasis.poly_vectors) == 20
    @test typeof(MyBasis.poly_functions[1]) <: Function
end    

@testset "PDEfind.jl" begin

    time_array, space_array = collect(0.0:0.03:0.4), collect(0.0:0.1:4.3)
    g = TimeSpaceGrid1D(time_array, space_array)
    data_array = time_array * space_array'
    ∂t = TimeDerivative(g, order=4)
    dt_data_array = ∂t(data_array) 
    dt_data_array_flat = reshape(dt_data_array, g.N) 
    ∂x = XDerivative(g; order=4) 
    max_derivative_degree = 3
    variables_matrix = GetVariablesMatrix(max_derivative_degree, data_array, g, ∂x)
    max_poly_degree = 2
    n_variables = max_derivative_degree + 1 
    MyBasis = PolynomialBasis(max_poly_degree, n_variables, skip_constant = false)
    Θ = evaluate(MyBasis, variables_matrix)

    xi, act_poly_vect = STRidge_cascade(Θ, dt_data_array_flat, MyBasis.poly_vectors, λ = 1e-5, tol = 0.0, iters = 4, verbose = false, normalize_columns = false)
    @test length(xi) == length(act_poly_vect)
    @test typeof(xi[1]) == Float64

    xi, act_poly_vect = TrainSTRidge(Θ, dt_data_array_flat, MyBasis.poly_vectors, λ = 1e-5, tol_multiplier=1.0, iters = 5, verbose = true, cond_number_multiplier=1e-6,  max_tol_iters = 100)
    @test length(xi) == length(act_poly_vect)
    
end