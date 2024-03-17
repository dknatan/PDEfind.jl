
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
    )

with `order` being the order of approximation of the finite difference method. 

# Fields 

* `M::T`

"""
struct TimeDerivative{T}
   M::T
end 

function TimeDerivative(
    g::TimeSpaceGrid1D; 
    order::Int = 2
)
    nt = g.nt
    dt = g.dt
    
    if order == 2
        M = diagm(1=>ones(nt-1)) - diagm(-1=>ones(nt-1))
        M[1,1] = -3.0
        M[1,2] = 4.0
        M[1,3] = -1.0   

        M[nt,nt] = 3.0
        M[nt,nt-1] = -4.0
        M[nt,nt-2] = 1.0
        return TimeDerivative(sparse(M./(2*dt)))
    elseif order == 4
        M = 8*(diagm(1=>ones(nt-1)) - diagm(-1=>ones(nt-1)))
        M+=-1*(diagm(2=>ones(nt-2)) - diagm(-2=>ones(nt-2)))
        M[1,1] = -25.0
        M[1,2] = 48.0
        M[1,3] = -36.0
        M[1,4] = 16.0
        M[1,5] = -3.0

        M[nt,nt] = 25.0
        M[nt,nt-1] = -48.0
        M[nt,nt-2] = 36.0
        M[nt,nt-3] = -16.0
        M[nt,nt-4] = 3.0

        M[2,1] = 0.0
        M[2,2] = -25.0
        M[2,3] = 48.0
        M[2,4] = -36.0
        M[2,5] = 16.0
        M[2,6] = -3.0

        M[nt-1,nt] = 0.0
        M[nt-1,nt-1] = 25.0
        M[nt-1,nt-2] = -48.0
        M[nt-1,nt-3] = 36.0
        M[nt-1,nt-4] = -16.0
        M[nt-1,nt-5] = 3.0
        return TimeDerivative(sparse(M./(12*dt)))
    else 
        throw("Only order 2 & 4 implemented")
    end 
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

with `order` being the order of approximation of the finite difference method and `boundary_condition` either "periodic" or "neumann" (no periodicity).

# Fields 

* `M::T`

"""
struct XDerivative{T}
    M::T
 end 
 
function XDerivative(g::TimeSpaceGrid1D; degree::Int = 1, order::Int = 2, boundary_condition::String = "neumann")
    nx = g.nx
    dx = g.dx
    
    if order == 2 && degree == 1
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
    elseif order == 2 && degree == 2
        M = diagm(1=>ones(nx-1)) + diagm(-1=>ones(nx-1)) - 2*diagm(0=>ones(nx))
        if boundary_condition == "periodic"
            M[1,nx] = -1.0
            M[nx,1] = 1.0
        elseif boundary_condition == "neumann"
            M[1,1] = 2.0
            M[1,2] = -5.0
            M[1,3] = 4.0   
            M[1,4] = -1.0   

            M[nx,nx] = 2.0
            M[nx,nx-1] = -5.0
            M[nx,nx-2] = 4.0
            M[nx,nx-3] = -1.0
        else
            throw("Only periodic & neumann BC implemented")
        end
        return XDerivative(sparse(M./dx^2))

    elseif order == 4 && degree == 1
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
