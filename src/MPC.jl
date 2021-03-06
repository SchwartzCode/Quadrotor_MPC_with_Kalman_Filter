# Taken from homework2-Q3

# All of the code in this cell is provided. There is no need to modify it.
"""
    MPCController

An MPC controller that uses a solver of type `S` to solve a QP at every iteration.

It will track the reference trajectory specified by `Xref`, `Uref` and `times` 
with an MPC horizon of `Nmpc`. It will track the terminal reference state if 
the horizon extends beyond the reference horizon.
"""
struct MPCController{S}
    P::SparseMatrixCSC{Float64,Int}
    q::Vector{Float64}
    A::SparseMatrixCSC{Float64,Int}
    lb::Vector{Float64}
    ub::Vector{Float64}
    Nmpc::Int
    solver::S
    Xref::Vector{Vector{Float64}}
    Uref::Vector{Vector{Float64}}
    times::Vector{Float64}
end


"""
    OSQPController(n,m,N,Nref,Nd)

Generate an `MPCController` that uses OSQP to solve the QP.
Initializes the controller with matrices consistent with `n` states,
`m` controls, and an MPC horizon of `N`, and `Nd` constraints. 

Use `Nref` to initialize a reference trajectory whose length may differ from the 
horizon length.
"""
function OSQPController(n::Integer, m::Integer, N::Integer, Nref::Integer=N, Nd::Integer=(N-1)*n)
    Np = (N-1)*(n+m)   # number of primals
    P = spzeros(Np,Np)
    q = zeros(Np)
    A = spzeros(Nd,Np)
    lb = zeros(Nd)
    ub = zeros(Nd)
    Xref = [zeros(n) for k = 1:Nref]
    Uref = [zeros(m) for k = 1:Nref]
    tref = zeros(Nref)
    solver = OSQP.Model()
    MPCController{OSQP.Model}(P,q, A,lb,ub, N, solver, Xref, Uref, tref)
end

isconstrained(ctrl::MPCController) = length(ctrl.lb) != (ctrl.Nmpc - 1) * length(ctrl.Xref[1])

"""
    buildQP!(ctrl, A,B,Q,R,Qf; kwargs...)

Build the QP matrices `P` and `A` for the MPC problem. Note that these matrices
should be constant between MPC iterations.

Any keyword arguments will be passed to `initialize_solver!`.
"""
function buildQP!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
    if isconstrained(ctrl)
        buildQP_constrained!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
    else
        buildQP_unconstrained!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
    end
end

"""
    updateQP!(ctrl::MPCController, x, time)

Update the vectors in the QP problem for the current state `x` and time `time`.
This should update `ctrl.q`, `ctrl.lb`, and `ctrl.ub`.
"""
function updateQP!(ctrl::MPCController, A,B, x, time)
    if isconstrained(ctrl)
        updateQP_constrained!(ctrl, A,B, x, time)
    else
        updateQP_unconstrained!(ctrl, x, time)
    end
end


"""
    initialize_solver!(ctrl::MPCController; kwargs...)

Initialize the internal solver once the QP matrices are initialized in the 
controller.
"""
function initialize_solver!(ctrl::MPCController{OSQP.Model}; tol=1e-6, verbose=false)
    OSQP.setup!(ctrl.solver, P=ctrl.P, q=ctrl.q, A=ctrl.A, l=ctrl.lb, u=ctrl.ub, 
        verbose=verbose, eps_rel=tol, eps_abs=tol, polish=1)
end


"""
Including build/updateQP_unconstrained in case control from code borrowed from homework2-Q3 calls
these functions
"""
function buildQP_unconstrained!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
    
    println("ERROR! buildQP_unconcstrained! not implemented")
    error("unimplemented")
    
    return nothing
end


function updateQP_unconstrained!(ctrl::MPCController, x, time)
    
    println("ERROR! updateQP_unconcstrained! not implemented")
    error("unimplemented")
        
    return nothing
end

"""
    get_k(ctrl, t)

Get the time index corresponding to time `t`. 
Useful for implementing zero-order hold control.
Uses binary search to find the time index.
"""
get_k(controller, t) = searchsortedlast(controller.times, t)

"""
    get_control(ctrl::MPCController, x, t)

Get the control from the MPC solver by solving the QP. 
If you want to use your own QP solver, you'll need to change this
method.
"""
function get_control(ctrl::MPCController{OSQP.Model}, A_traj, B_traj, x, time; wind_correction=true)
    
    # Update the QP
    updateQP!(ctrl, A_traj, B_traj, x, time)
    OSQP.update!(ctrl.solver, q=ctrl.q, l=ctrl.lb, u=ctrl.ub)
    
    # Solve QP
    results = OSQP.solve!(ctrl.solver)
    ??u = results.x[1:4]
    
    k = get_k(ctrl, time)
    umax = 15.0
    umin = -2.0
    
    u_corrected_pre_wind = ctrl.Uref[k] + ??u
    clamp!(u_corrected_pre_wind, umin, umax)
    
    A,B = KF_dynamics_jacobians(x, u_corrected_pre_wind, dt)
    du_wind = get_wind_correction(x, B, dt) #note: this function is in dynamics.jl
    
    u_corrected = ctrl.Uref[k] + ??u
    
    if wind_correction
        u_corrected += du_wind
    end
    
    clamp!(u_corrected, umin, umax)
    
    println("wind correction: ", u_corrected - u_corrected_pre_wind)
    
    return u_corrected
    
end


"""
    buildQP!(ctrl, A,B,Q,R,Qf; kwargs...)

Build the QP matrices `P` and `A` for the MPC problem. Note that these matrices
should be constant between MPC iterations.

Any keyword arguments will be passed to `initialize_solver!`.
"""
function buildQP_constrained!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
    # TODO: Implement this method to build the QP matrices
    
    x_size = size(ctrl.Xref[1])[1] - 1 # -1 to account for quaternion -> 3-param attitude
    u_size = size(ctrl.Uref[1])[1]
    
    # useful variables
    N = ctrl.Nmpc
    R_size = size(R)
    Q_size = size(Q)
    P_rows = (N-1)*(u_size + x_size)
    P_cols = (N-1)*(u_size + x_size)
    
    # populate first row of A matrix
    A_constraint = [B[1] -I(x_size) zeros(x_size, P_cols - u_size - x_size)]
    
    # set first row of P
    P = [R zeros(u_size, P_cols - u_size)]
    curr_pos = u_size
    
#     Inequality constraint matrices that will end up in A
    # TODO: replace all this matrix stuff with kron calls
    u1_selector = zeros(u_size)
    u1_selector[1] = 1
    A_u1_constraint = [u1_selector' zeros(P_cols - curr_pos)']
    
    u2_selector = zeros(u_size)
    u2_selector[2] = 1
    A_u2_constraint = [u2_selector' zeros(P_cols - curr_pos)']
    
    u3_selector = zeros(u_size)
    u3_selector[3] = 1
    A_u3_constraint = [u3_selector' zeros(P_cols - curr_pos)']
    
    u4_selector = zeros(u_size)
    u4_selector[4] = 1
    A_u4_constraint = [u4_selector' zeros(P_cols - curr_pos)']
        
    for i=2:N-1
        # Rows of A constraint matrix
        new_A_const_row = [zeros(x_size, curr_pos) A[i] B[i] -I(x_size) zeros(x_size, P_cols - curr_pos - 2*x_size - u_size)]
        A_constraint = [A_constraint; new_A_const_row]
        
        # ==== P rows ====
        # Q row
        new_P_row = [zeros(x_size, curr_pos) Q zeros(x_size, P_cols - (curr_pos + x_size))]
        curr_pos += x_size
        P = [P; new_P_row]
            
        # R row
        new_P_row = [zeros(u_size, curr_pos) R zeros(u_size, P_cols - (curr_pos + u_size))]
        curr_pos += u_size
        P = [P; new_P_row]
        
        A_u1_constraint = [A_u1_constraint; zeros(curr_pos - u_size)' u1_selector' zeros(P_cols - curr_pos)']
        A_u2_constraint = [A_u2_constraint; zeros(curr_pos - u_size)' u2_selector' zeros(P_cols - curr_pos)']
        A_u3_constraint = [A_u3_constraint; zeros(curr_pos - u_size)' u3_selector' zeros(P_cols - curr_pos)']
        A_u4_constraint = [A_u4_constraint; zeros(curr_pos - u_size)' u4_selector' zeros(P_cols - curr_pos)']
    end
    
    final_P_row = [zeros(x_size, P_cols - x_size) Qf]
    P = [P; final_P_row]
    
    ctrl.P .= P
    
    A_constraint = [A_constraint; A_u1_constraint; A_u2_constraint; A_u3_constraint; A_u4_constraint]
    ctrl.A .= A_constraint    
    
    # Initialize the included solver
    #    If you want to use your QP solver, you should write your own
    #    method for this function
    initialize_solver!(ctrl; kwargs...)
    return nothing
end

"""
    update_QP!(ctrl::MPCController, x, time)

Update the vectors in the QP problem for the current state `x` and time `time`.
This should update `ctrl.q`, `ctrl.lb`, and `ctrl.ub`.
"""
function updateQP_constrained!(ctrl::MPCController, A, B, x, time)
    k = get_k(ctrl, time)
    x_size = size(ctrl.Xref[1])[1] - 1
    u_size = size(ctrl.Uref[1])[1]
    N = ctrl.Nmpc
    
    # TODO: define these not in place
    xeq = ctrl.Xref[end]
    ueq = ctrl.Uref[1]
    
    R = ctrl.P[begin:u_size, begin:u_size]
    Q = ctrl.P[u_size+1:u_size+x_size, u_size+1:u_size+x_size]
    Qf = ctrl.P[1+end-x_size:end, 1+end-x_size:end]
        
    state_bounds = zeros(x_size*(ctrl.Nmpc-1))
    d?? = ??(quat_L(ctrl.Xref[k][4:7])' * x[4:7])
    dX = vcat([(x - ctrl.Xref[k])[1:3]' d??' (x - ctrl.Xref[k])[8:end]'])'
            
#     A_o, B_o = dynamics_jacobians(x,Uref[k],dt)
#     J_attitude = attitude_jacobian(x)
#     A_o = J_attitude'*A_o*J_attitude
#     B_o = J_attitude'*B_o
    
    state_bounds[begin:size(A[k])[1]] = -A[k]*dX
    
    # TODO: update along Uref if not all uhover
    thrust_ub = (15.0-Uref[1][1]) * ones(ctrl.Nmpc-1)
    thrust_lb = (-2.0-Uref[1][1]) * ones(ctrl.Nmpc-1)
#     thrust_ub = 50000 * ones(ctrl.Nmpc-1)
#     thrust_lb = -20000 * ones(ctrl.Nmpc-1)
    
    ub = [state_bounds; thrust_ub; thrust_ub; thrust_ub; thrust_ub]
    lb = [state_bounds; thrust_lb; thrust_lb; thrust_lb; thrust_lb]
    
    ctrl.ub .= ub
    ctrl.lb .= lb
    
    P_rows = (N-1)*(u_size + x_size)
    P_cols = (N-1)*(u_size + x_size)
                
    # Update jacobians in ctrl.A
    ctrl.A[1:x_size, :] .= [B[k] -I(x_size) zeros(x_size, P_cols - u_size - x_size)]
    
    curr_pos = u_size
    
    for i=1:N-2
        # Rows of A constraint matrix
        if k+i > length(A)
            new_A_const_row = zeros(x_size, P_cols)
        else
            new_A_const_row = [zeros(x_size, curr_pos) A[k+i] B[k+i] -I(x_size) zeros(x_size, P_cols - curr_pos - 2*x_size - u_size)]
        end
        ctrl.A[(1:x_size) .+ (i*x_size),:] .= new_A_const_row

        curr_pos += x_size + u_size 
    end
            
                
    return nothing
end


"""
Create MPC object and solve
"""
function build_MPC_QP(Xref, Uref, tref, A, B, Q, R, Qf)
    # Initialize the constrained MPC controller
    n = size(B[1])[1]
    m = size(B[1])[2]

    Nd = (Nmpc-1)*(n+m)
    mpc = OSQPController(n, m, Nmpc, length(Xref), Nd)
    mpc.Xref .= Xref
    mpc.Uref .= Uref
    mpc.times .= tref
    
    # TODO: do we need/want to adjust A and B along the trajectory?
    buildQP!(mpc, A, B, Q, R, Qf, tol=1e-2, verbose=false)
    
    return mpc
end