using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets, Plots

function GetExampleData(
    # pdetype::String,
    # domain,
    # initial_condition,
    # boundary_condition,
    # parameters,
    # noise
)    
    
    # Parameters, variables, and derivatives
    @parameters t x
    @variables u(..)
    Dt = Differential(t)
    Dxx = Differential(x)^2

    # 1D PDE and boundary conditions
    eq  = Dt(u(t, x)) ~ Dxx(u(t, x))
    bcs = [u(0, x) ~ cos(x),
            u(t, 0) ~ exp(-t),
            u(t, 1) ~ exp(-t) * cos(1)]

    # Space and time domains
    domains = [t ∈ Interval(0.0, 1.0),
            x ∈ Interval(0.0, 1.0)]

    # PDE system
    @named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])

    # Method of lines discretization
    dx = 0.1
    order = 2 # TODO why is this needed
    discretization = MOLFiniteDifference([x => dx], t)

    # Convert the PDE problem into an ODE problem
    prob = discretize(pdesys,discretization)

    # Solve ODE problem
    sol = solve(prob, Tsit5(), saveat=0.2)
    

    # Method of Manufactured Solutions: exact solution
    u_exact = (x,t) -> exp.(-t) * cos.(x)

    # Plot results and compare with exact solution
    discrete_x = sol[x]
    discrete_t = sol[t]
    solu = sol[u(t, x)]

    plt = plot()

    for i in eachindex(discrete_t)
        plot!(discrete_x, solu[i, :], label="Numerical, t=$(discrete_t[i])")
        scatter!(discrete_x, u_exact(discrete_x, discrete_t[i]), label="Exact, t=$(discrete_t[i])")
    end
    plt

    return "hi"
end
