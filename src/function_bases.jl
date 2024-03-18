abstract type AbstractBasis end

"""
    PolynomialBasis <: AbstractBasis

Basis of polynomial functions up to degree `max_degree` of `n_variables` variables.

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
    skip_constant::Bool=false
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

Construct all polynomials of degree exactly `degree` of `n_variables` variables.
Returns array of arrays of length `n_variables`, `i`-th entry representing the power of `i`-th variable. 
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
                push!(poly_collection, p_candidate) 
            end
        end
        return poly_collection
    end
end

"""
    all_polynomial_vectors(max_degree::Int, n_variables::Int)

Construct all polynomials up to degree `max_degree` of `n_variables` variables. 
Returns array of arrays of length `n_variables`, `i`-th entry representing the power of `i`-th variable. 
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

Construct all polynomials up to degree `max_degree` of `n_variables` variables. 
Returns array of functions of `n_variables` variables. Only works for 5 variables or less.
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
        throw("More than 5 variables not implemented. ")
    end
end

"""
    function evaluate(basis::PolynomialBasis, variables_array)

Evaluate basis functions at different values of variables.
If `basis` includes M functions (of n variables) and the `variables_array` is an N × n array, the output is a N × M array.

"""
function evaluate(basis::PolynomialBasis, variables_array)
    return hcat([polynomial.(variables_array...) for polynomial in basis.poly_functions]...)
end
