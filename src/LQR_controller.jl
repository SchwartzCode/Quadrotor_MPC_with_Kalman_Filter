# from homework 3 problem 1
# I think all these funcs are for iLQR... probably don't need them?

function stage_cost(x,u,xref,uref)
    # LQR cost at each knot point (depends on both x and u)
    x_diff = x - xref
    u_diff = u - uref
    
    J = 0.5 * x_diff' * Q * x_diff + 0.5 * u_diff' * R * u_diff 
    
    return J
end

function term_cost(x,xref)
    # LQR terminal cost (depends on just x)
    x_diff = x - xref
    
    J = 0.5 * x_diff' * Qf * x_diff
    
    return J
end

function trajectory_cost(X,U,Xref,Uref)
    # calculate the cost of a given trajectory 
    J = 0.0
    for i = 1:size(U)[1]
        J += stage_cost(X[i], U[i], Xref[i], Uref[i])
    end
    J += term_cost(X[end], Xref[end])
    
    return J
end
        
function stage_cost_expansion(x,u,xref,uref)
    # if the stage cost function is J, return the following derivatives:
    # ∇²ₓJ,  ∇ₓJ, ∇²ᵤJ, ∇ᵤJ
    Jxx = Q #zeros(nx,nx)
    Jx = Q*(x - xref) #zeros(nx)
    Juu = R #zeros(nu,nu)
    Ju = R * (u - uref) #zeros(nu)
    
    return Jxx, Jx, Juu, Ju
end

function term_cost_expansion(x,xref)
    # if the terminal cost function is J, return the following derivatives:
    # ∇²ₓJ,  ∇ₓJ
    Jxx = Qf #zeros(nx,nx)
    Jx = Qf*(x - xref) #zeros(nx)
    
    return Jxx, Jx
end

# From homework3 problem 2
function lqr(A,B,Q,R; max_iters=200, tol=1e-6)
    
    n,m = size(B)
    K = zeros(m,n)
    K_prev = copy(K)
    P = copy(Q)
    
    for k = 1:max_iters        
        K .= (R + B'P*B) \ (B'P*A)
        P .= Q + A'P*A - A'P*B*K
        
#         println("\t", k, ":  ", norm(K-K_prev,Inf))
        if norm(K-K_prev,Inf) < tol
            println("LQR converged in $k iterations")
            return K
        end
        K_prev .= K
    end
    println("no LQR convergence, final vals:\t", K)
    return K * NaN
end