module PDEfind

using SparseArrays, LinearAlgebra 

# export GetExampleData
# include("../examples/PDEexamples.jl")
# # Write your package code here.

export TimeSpaceGrid1D, TimeDerivative, XDerivative, GetVariablesMatrix, AbstractBasis, PolynomialBasis, evaluate, STRidge_cascade, TrainSTRidge, lambda_sweep

include("grid_operations.jl")
include("function_bases.jl")


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
    poly_vectors; 
    λ::Float64 = 0.0, 
    tol::Float64 = 0.01, 
    iters::Int64 = 1,
    verbose::Bool = false,
    normalize_columns::Bool = true
    )
    
    function _print_tool(active_poly_vectors, ξ, tol, iters)
        println("Iteration: $iters, threshold: $tol. ")

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
            push!(parameter_names, "\\xi_{$(poly_vector...)}")
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

    function STRidge_recursive(theta::Matrix{Float64}, dt_data_array_flat::Vector{Float64}, active_poly_vectors::Array, norm_facts::Vector{Float64}; λ::Float64 = 0.0, tol::Float64 = 0.02, iters::Int64 = 1)
        # Sanity checks
        @assert size(theta, 1) == length(dt_data_array_flat)
        @assert size(theta, 2) == length(active_poly_vectors)
        @assert size(theta, 2) == length(norm_facts)

        # Ridge regression
        N = size(theta, 1)
        # ξ = inv(theta' * theta + λ * diagm(0=>ones(size(theta, 2)))) * theta' * dt_data_array_flat
        ξ_norm = (theta' * theta ./ N + λ .* diagm(0=>ones(size(theta, 2)))) \ theta' * dt_data_array_flat ./ N
        ξ = ξ_norm ./ norm_facts
        # Threshold
        bigcoeffs = abs.(ξ) .> tol
        
        if verbose
            _print_tool(active_poly_vectors, ξ, tol, iters)
        end 

        # Base case
        if iters == 1
            return ξ[bigcoeffs], active_poly_vectors[bigcoeffs]

        # Recursive call
        else
            return STRidge_recursive(theta[:, bigcoeffs], dt_data_array_flat, active_poly_vectors[bigcoeffs], norm_facts[bigcoeffs]; λ=λ, tol=tol, iters=iters-1)
        end
    end

    if normalize_columns
        N = size(theta, 1)
        means = sum(theta, dims = 1) ./ N
        vars = sum(theta .^ 2 .- means .^ 2, dims = 1) ./ N
        norm_facts = sqrt.(vars)
        if vars[1] == 0.0 # fix var(1, 1, ...) = 0 problem
            norm_facts[1] = 1.0
        end

        theta_norm = theta ./ norm_facts 
        
        norm_facts = reshape(norm_facts, size(norm_facts, 2))
    else
        theta_norm = copy(theta)
        norm_facts = ones(size(theta, 2))
    end 


    # theta_final, active_poly_vectors_final =  STRidge_recursive(theta, dt_data_array_flat, poly_vectors; λ=λ, tol=tol, iters=iters)
    
    # # Final non-ridge regression
    # ξ = inv(theta_final' * theta_final) * theta_final' * dt_data_array_flat
    
    return STRidge_recursive(theta_norm, dt_data_array_flat, poly_vectors, norm_facts; λ=λ, tol=tol, iters=iters)
end


""" 
    TrainSTRidge(
        theta::Matrix{Float64}, 
        dt_data_array_flat::Vector{Float64}, 
        poly_vectors; 
        λ::Float64 = 1e-4, 
        tol_multiplier::Float64 = 1e2,
        iters::Int64 = 1,
        verbose::Bool = true,
        cond_number_multiplier::Float64 = 1e-2,
        max_tol_iters::Int64 = 100
    )
    Conduct sequential thresholded ridge regression. See algorithm 2 in Rudy et al and the .pdf file attached by the author.
    Returns optimal ξ. 
"""
function TrainSTRidge(
    theta::Matrix{Float64}, 
    dt_data_array_flat::Vector{Float64}, 
    poly_vectors; 
    λ::Float64 = 1e-4, 
    tol_multiplier::Float64 = 1e2,
    iters::Int64 = 1,
    verbose::Bool = true,
    cond_number_multiplier::Float64 = 1e-2,
    max_tol_iters::Int64 = 100
)
    @assert size(theta, 1) == length(dt_data_array_flat)
    @assert size(theta, 2) == length(poly_vectors)

    # split theta, dt_data_array_flat to train&test subarrays
    N = size(theta, 1)
    N_train = Int64(round(0.8 * N))
    N_test = N - N_train
    theta_train, theta_test = theta[1:N_train, :], theta[N_train:end, :]
    dt_train, dt_test = dt_data_array_flat[1:N_train], dt_data_array_flat[N_train:end]

    # empirical l0 penalty
    cond_number = cond(theta, 2) # not 100% sure this condition number makes sense
    η = cond_number_multiplier * cond_number
    
    if verbose
        println("Split input data to $N_train train and $N_test test samples. ")
        println("Using η = $cond_number_multiplier * $cond_number. ")
    end

    # baseline predictor
    ξ_base, active_poly_vectors_base = STRidge_cascade(theta_train, dt_train, poly_vectors, λ=λ, tol=0.0, iters=1, verbose=false)
    tol = sum(abs.(ξ_base)) + 1.0 
    ξ_best, active_poly_vectors_best = STRidge_cascade(theta_train, dt_train, poly_vectors, λ=λ, tol=tol, iters=2, verbose=false)
    error_best = norm(dt_test, 2)^2
    @assert length(ξ_best) == 0

    # search through Threshold values
    for tol_iter in 1:max_tol_iters
        # take step
        tol /= 1. + tol_multiplier
        
        # get xi from train set with current iters and evaluate error
        ξ, active_poly_vectors = STRidge_cascade(theta_train, dt_train, poly_vectors, λ=λ, tol=tol, iters=iters, verbose=false)
        chosen_mask = poly_vectors .∈ Ref(active_poly_vectors)
        error = norm(theta_test[:, chosen_mask] * ξ - dt_test, 2)^2 + η * norm(ξ, 0)^2

        if length(ξ) - length(ξ_best) > 1
            # step back and decrease stepsize multiplier
            tol *= 1. + tol_multiplier
            tol_multiplier /= 2.
            if verbose
                println("tol_iter = $tol_iter: Stepsize too big. Decreased tol_multiplier to $tol_multiplier. ")
            end
         
        else
            if error > error_best
                # done
                if verbose
                    println("tol_iter = $tol_iter: Found optimal fit with polynomials $active_poly_vectors_best, xi = $ξ_best")
                end
                break
            else
                # take step
                ξ_best, active_poly_vectors_best, error_best = ξ, active_poly_vectors, error
                # tol_multiplier *= 2. # just a speed up i think
                if verbose 
                    println("tol_iter = $tol_iter: Found better/equal fit with polynomials $active_poly_vectors_best, xi = $ξ_best. Decreased Threshold to $tol. ")
                end    
            end
        end
    end
    return ξ_best, active_poly_vectors_best
end

"""
    function lambda_sweep(
        lambda_range,
        Θ::Matrix{Float64}, 
        dt_data_array_flat::Vector{Float64},
        poly_vectors;
        iters::Int64 = 1,
        cond_number_multiplier::Float64 = 1e-6
    )

    Evaluate optimal ξ at all values of `lambda_range`. Returns array of size `(length(poly_vectors), 2)` with optimal (λ, ξ_i) values for each possible polynomial. 
"""
function lambda_sweep(
    lambda_range,
    Θ::Matrix{Float64}, 
    dt_data_array_flat::Vector{Float64},
    poly_vectors;
    iters::Int64 = 1,
    cond_number_multiplier::Float64 = 1e-6
)
    res = Array{Any}(undef, length(poly_vectors), 2)
    for λ ∈ lambda_range
        ξ_array, active_poly_vectors = TrainSTRidge(Θ, dt_data_array_flat, poly_vectors, λ=λ, tol_multiplier=1.0, iters=iters, verbose=false, cond_number_multiplier=cond_number_multiplier)
        for (ξ,poly_vector) ∈ zip(ξ_array, active_poly_vectors) 
            poly_index = findall(x -> x == poly_vector, poly_vectors)[1]
            try
                push!(res[poly_index, 1], λ)
                push!(res[poly_index, 2], ξ)    
            catch
                res[poly_index, 1] = [λ]
                res[poly_index, 2] = [ξ]
            end
        end
    end
    return res
end

end
