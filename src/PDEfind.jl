module PDEfind

using SparseArrays, LinearAlgebra 

# export GetExampleData
# include("../examples/PDEexamples.jl")
# # Write your package code here.

export TimeSpaceGrid1D, TimeDerivative, XDerivative, GetVariablesMatrix, AbstractBasis, PolynomialBasis, evaluate, STRidge_cascade

"""
    TimeSpaceGrid1D{T}

Dicretization grid for time + 1D space finite difference schemes. 

# Initialization 

    TimeSpaceGrid1D(time_array, space_array)

with both assumed to have constant spacing 

# Fields 

* `t::AbstractArray{T}`
* `x::AbstractArray{T}`
* `dt::T`
* `dx::T`
* `nt::Int`
* `nx::Int`
* `N::Int`
* `t_grid::AbstractMatrix{T}`
* `x_grid::AbstractMatrix{T}`
"""
struct TimeSpaceGrid1D{T}
    t::AbstractArray{T} # this is a n_t long range
    x::AbstractArray{T}
    dt::T
    dx::T
    nt::Int # number of time grid points
    nx::Int # number of space grid points
    N::Int # total number of grid points 
    t_grid::AbstractMatrix{T} # this is a nt×nx matrix
    x_grid::AbstractMatrix{T} # this is a nt×nx matrix
end 

function TimeSpaceGrid1D(t::AbstractArray{T}, x::AbstractArray{T}) where T

    dt = abs(t[2] - t[1])
    dx = abs(x[2] - x[1])

    nt = length(t)
    nx = length(x)
    
    N = nt * nx
    t_grid = ones(nx) * t' # meshgrid 
    x_grid = x * ones(nt)' # meshgrid
    
    return TimeSpaceGrid1D(t, x, dt, dx, nt, nx, N, t_grid, x_grid)
end 


"""
    TimeDerivative{T}

Finite difference operator for time derivative on a grid. 

# Initialization 

    TimeDerivative(
        g::TimeSpaceGrid1D; 
        order::Int = 2,
        boundary_condition::String = "neumann"
    )

with ... 

# Fields 

* `M::T`

# Usage
    TODO
"""
struct TimeDerivative{T}
   M::T
end 

function TimeDerivative(
    g::TimeSpaceGrid1D; 
    order::Int = 2,
    boundary_condition::String = "neumann"
)
    nt = g.nt
    dt = g.dt
    
    if order == 2
        M = diagm(1=>ones(nt-1)) - diagm(-1=>ones(nt-1))
        if boundary_condition == "neumann"
            M[1,1] = -2.0
            M[1,2] = 2.0
            M[nt,nt-1] = -2.0
            M[nt,nt] = 2.0
        else 
            throw("only neumann BC (time)")
        end
        TimeDerivative(sparse(M./(2*dt)))
    else 
        throw("Only order 2 implemented")
    end 
    return TimeDerivative(M)
end

function (∂t::TimeDerivative)(data::AbstractMatrix)
    return ∂t.M * data 
end 


"""
    XDerivative{T}

Finite difference operator for (1st) spatial derivative on a grid. 

# Initialization 

    XDerivative(
        g::TimeSpaceGrid1D; 
        order::Int = 2,
        boundary_condition::String = "neumann"
    )

with ... 

# Fields 

* `M::T`

# Usage
    TODO
"""
struct XDerivative{T}
    M::T
 end 
 
function XDerivative(g::TimeSpaceGrid1D; order::Int = 2, boundary_condition::String = "neumann")
    nx = g.nx
    dx = g.dx
    
    if order == 2
        M = diagm(1=>ones(nx-1)) - diagm(-1=>ones(nx-1))
        if boundary_condition == "periodic"
            M[1,nx] = -1.0
            M[nx,1] = 1.0
        elseif boundary_condition == "neumann"
            M[1,1] = -3.0
            M[1,2] = 4.0
            M[1,3] = -1.0   

            M[nx,nx] = 3.0
            M[nx,nx-1] = -4.0
            M[nx,nx-2] = 1.0
        else
            throw("Only periodic & neumann BC implemented")
        end
        return XDerivative(sparse(M./(2*dx)))
    elseif order == 4
        M = 8*(diagm(1=>ones(nx-1)) - diagm(-1=>ones(nx-1)))
        M+=-1*(diagm(2=>ones(nx-2)) - diagm(-2=>ones(nx-2)))
        if boundary_condition == "periodic"
            M[1,nx] = -8.0
            M[nx,1] = 8.0
            M[2,nx] = 1.0
            M[nx,2] = -1.0
        elseif boundary_condition == "neumann"
            M[1,1] = -25.0
            M[1,2] = 48.0
            M[1,3] = -36.0
            M[1,4] = 16.0
            M[1,5] = -3.0

            M[nx,nx] = 25.0
            M[nx,nx-1] = -48.0
            M[nx,nx-2] = 36.0
            M[nx,nx-3] = -16.0
            M[nx,nx-4] = 3.0

            M[2,1] = 0.0
            M[2,2] = -25.0
            M[2,3] = 48.0
            M[2,4] = -36.0
            M[2,5] = 16.0
            M[2,6] = -3.0

            M[nx-1,nx] = 0.0
            M[nx-1,nx-1] = 25.0
            M[nx-1,nx-2] = -48.0
            M[nx-1,nx-3] = 36.0
            M[nx-1,nx-4] = -16.0
            M[nx-1,nx-5] = 3.0
        else
            throw("Only periodic & neumann BC implemented")
        end
        
        return XDerivative(sparse(M./(12*dx)))
    else 
        throw("Only order 2 & 4 implemented")
    end
end

function (∂x::XDerivative)(data::AbstractMatrix)
    return (data * ∂x.M')
end 

"""
    GetVariablesMatrix(
        max_derivatives::Int,
        data_array::Matrix,
        g::TimeSpaceGrid1D,
        ∂x::XDerivative
)   

Create a matrix (nt*nx, max_derivatives+1) containing all j-derivatives evaluated at all time-space coordinates i.
Note: This is NOT the Θ matrix. 
"""
function GetVariablesMatrix(
    max_derivatives::Int,
    data_array::Matrix,
    g::TimeSpaceGrid1D,
    ∂x::XDerivative
)   
    x_derivatives_flat = [reshape(data_array, g.N)] 
    current_array = data_array
    for i in 1:max_derivatives
        current_array = ∂x(current_array)
        push!(x_derivatives_flat, reshape(current_array, g.N))
    end
    return x_derivatives_flat
end


abstract type AbstractBasis end

"""
    PolynomialBasis <: AbstractBasis

Basis of polynomial functions up to max_degree of n_variables 

# Initialization 

    PolynomialBasis(
        max_degree::Int,
        n_variables::Int
    )

with ... 

# Fields 

* `max_degree::Int`
* `n_variables::Int`
* `poly_vectors::Array`
* `poly_functions::Array{Function}`

# Usage
    TODO
"""
struct PolynomialBasis <: AbstractBasis
    max_degree::Int
    n_variables::Int
    poly_vectors::Array
    poly_functions::Array{Function}
end

function PolynomialBasis(
    max_degree::Int,
    n_variables::Int;
    skip_constant::Bool=true
)
    poly_vectors = all_polynomial_vectors(max_degree, n_variables)
    poly_functions = all_polynomial_functions(max_degree, n_variables)
    if skip_constant
        poly_vectors = poly_vectors[2:end]
        poly_functions = poly_functions[2:end]
    end
    return PolynomialBasis(max_degree, n_variables, poly_vectors, poly_functions)
end

"""
    polynomial_vectors(degree::Int, n_variables::Int)

Returns all polynomials of degree exactly degree of n_variables.
Returns array of arrays of length n_variables, i-th entry representing the power of i-th variable. 
"""
function polynomial_vectors(degree::Int, n_variables::Int)
    if degree == 0
        return Any[zeros(Int64, n_variables)]
    elseif degree == 1
        poly_collection_deg1 = []
        for i in 1:n_variables
            poly_1 = zeros(Int64, n_variables)
            poly_1[i] = 1
            push!(poly_collection_deg1, poly_1)
        end
        return poly_collection_deg1
    else
        poly_collection_prev = polynomial_vectors(degree-1, n_variables)
        poly_collection_deg1 = polynomial_vectors(1, n_variables)
        poly_collection = []
        for poly_prev in poly_collection_prev, poly_1 in poly_collection_deg1
            p_candidate = poly_prev + poly_1
            if p_candidate in poly_collection
                continue
            else
                push!(poly_collection,p_candidate) 
            end
        end
        return poly_collection
    end
end

"""
    all_polynomial_vectors(max_degree::Int, n_variables::Int)

Returns all polynomials up to degree max_degree. 
Returns array of arrays of length n_variables, i-th entry representing the power of i-th variable. 
"""
function all_polynomial_vectors(max_degree::Int, n_variables::Int)
    all_poly_collection = []
    for i in 0:max_degree
        for poly_collection in polynomial_vectors(i,n_variables)
            push!(all_poly_collection, poly_collection)
        end
    end

    return all_poly_collection
end

"""
    all_polynomial_functions(max_degree::Int, n_variables::Int)

Returns all polynomials up to degree max_degree. 
Returns array of functions of n_variables. Only works for 5 variables or less.
"""
function all_polynomial_functions(max_degree::Int, n_variables::Int)
    all_p = all_polynomial_vectors(max_degree, n_variables)
    if n_variables == 0
        throw("No variables no work!")
    elseif n_variables == 1    
        return [((x ) -> x^p[1]) for p in all_p]
    elseif n_variables == 2
        return [((x, y) -> x^p[1] * y^p[2]) for p in all_p]
    elseif n_variables == 3
        return [((x, y, z) -> x^p[1] * y^p[2] * z^p[3]) for p in all_p]
    elseif n_variables == 4
        return [((x, y, z, w) -> x^p[1] * y^p[2] * z^p[3] * w^p[4]) for p in all_p]
    elseif n_variables == 5
        return [((x, y, z, w, u) -> x^p[1] * y^p[2] * z^p[3] * w^p[4] * u^p[5]) for p in all_p]
    else
        throw("You have more than 5 variables. You must be crazy")
    end
end


function evaluate(basis::PolynomialBasis, variables_array; normalize_columns::Bool = false)
    theta = hcat([polynomial.(variables_array...) for polynomial in basis.poly_functions]...)
    if normalize_columns
        N = size(theta, 1)
        means = sum(theta, dims = 1) ./ N
        vars = sum(theta .^ 2 .- means .^ 2, dims = 1) ./ N
        norm_facts = sqrt.(vars)
        if vars[1] == 0.0 # fix var(1, 1, ...) = 0 problem
            norm_facts[1] = 1.0
        end
        theta ./= norm_facts    
        return theta, norm_facts

    else
        return theta, nothing
    end 
end



"""
    STRidge_cascade(
        theta::Matrix{Float64}, 
        dt_data_array_flat::Vector{Float64}, 
        MyBasis::PolynomialBasis; 
        λ::Float64 = 0.0, 
        tol::Float64 = 0.01, 
        iters::Int64 = 1,
        verbose = false
    )    
    Conduct sequential thresholded ridge regression. See algorithm 1 in Rudy et al.
    Returns optimal ξ. 
"""
function STRidge_cascade(
    theta::Matrix{Float64}, 
    dt_data_array_flat::Vector{Float64}, 
    MyBasis::PolynomialBasis; 
    λ::Float64 = 0.0, 
    tol::Float64 = 0.01, 
    iters::Int64 = 1,
    verbose = false
)
    
    function _print_tool(active_poly_vectors, ξ, tol)
        # infer number of variables
        n_variables = length(active_poly_vectors[1])

        # construct all variable names
        variable_names = ["u"]
        for i in 1:(n_variables-1)
            push!(variable_names, "\\partial_{" * "x"^i * "} u")
        end

        # construct all (properly colored) parameter names
        parameter_names = []
        for (i,poly_vector) in enumerate(active_poly_vectors)
            push!(parameter_names, "p_{$(poly_vector...)}")
        end

        # rounding the coefficients
        decimal_places = length(split(string(tol), ".")[2]) + 1

        # initialize output string(s)
        output_string = "Current fit: \n \$ \\\\ \\partial_t u = "
        num_vals = ""
        bigcoeffs = abs.(ξ) .> tol

        # loop over active terms
        for (i,poly_vector) in enumerate(active_poly_vectors)
            term_to_add = ""
            # add parameter
            term_to_add *= parameter_names[i]
            # loop over factors of term
            for (j,poly_power) in enumerate(poly_vector) 
                # if clause to take care of parenthesis
                if poly_power == 1
                    term_to_add *= variable_names[j]
                elseif poly_power > 1
                    term_to_add *= "(" * variable_names[j] * ")^" * string(poly_power)
                end
            end
            # decide color
            term_to_add = bigcoeffs[i] ? term_to_add : "\\textcolor{grey}{" * term_to_add * "}"
            output_string *= term_to_add
            output_string *= " + "


            # parameter value (properly rounded)
            num_val_to_add = parameter_names[i] * " = $(round(ξ[i],digits=decimal_places)) \\\\ "
            # decide color
            num_vals *= bigcoeffs[i] ? num_val_to_add : "\\textcolor{grey}{" * num_val_to_add * "}"
        end

        # remove the final "+" sign
        if output_string[end-2:end] == " + "
            output_string = output_string[1:end-2]
        end

        # append numerical values
        output_string *= "\\\\ \$ with: \$ \\\\" * num_vals * "\$ \n - - -"

        display("text/markdown", output_string)
    end

    function STRidge_recursive(theta::Matrix{Float64}, dt_data_array_flat::Vector{Float64}, active_poly_vectors::Array; λ::Float64 = 0.0, tol::Float64 = 0.02, iters::Int64 = 1)
        # Sanity checks
        @assert size(theta, 1) == length(dt_data_array_flat)
        @assert size(theta, 2) == length(active_poly_vectors)

        # Ridge regression
        N = size(theta, 1)
        ξ = inv(theta' * theta ./ N + λ * diagm(0=>ones(size(theta, 2)))) * theta' * dt_data_array_flat ./ N
        
        # Threshold
        bigcoeffs = abs.(ξ) .> tol
        
        if verbose
            println("Iterations to go: $iters, Tolerance: $tol, λ = $λ")
            _print_tool(active_poly_vectors, ξ, tol)
        end

        # Base case
        if iters == 1
            return ξ # theta[:, bigcoeffs], active_poly_vectors[bigcoeffs]

        # Recursive call
        else
            return STRidge_recursive(theta[:, bigcoeffs], dt_data_array_flat, active_poly_vectors[bigcoeffs]; λ=λ, tol=tol, iters=iters-1)
        end
    end

    # theta_final, active_poly_vectors_final =  STRidge_recursive(theta, dt_data_array_flat, MyBasis.poly_vectors; λ=λ, tol=tol, iters=iters)
    
    # # Final non-ridge regression
    # ξ = inv(theta_final' * theta_final) * theta_final' * dt_data_array_flat
    
    return STRidge_recursive(theta, dt_data_array_flat, MyBasis.poly_vectors; λ=λ, tol=tol, iters=iters)
end


# function remove_polynomial(MyBasis::PolynomialBasis, poly_vector)
#     if poly_vector in MyBasis.poly_vectors
#         index = findall(x->x==poly_vector, MyBasis.poly_vectors)
#         MyBasis.poly_vectors = deleteat!(MyBasis.poly_vectors, index)    
#         MyBasis.poly_functions = deleteat!(MyBasis.poly_functions, index)    
#     else
#         println("The polynomial vector $poly_vector not found in current basis.")
#     end 
# end

# function remove_polynomials(MyBasis::PolynomialBasis, poly_vector_array)
#     for poly_vector in poly_vector_array
#         remove_polynomial(MyBasis, poly_vector)
#     end 
# end

end
