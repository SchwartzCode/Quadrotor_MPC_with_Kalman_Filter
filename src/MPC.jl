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
function updateQP!(ctrl::MPCController, x, time)
    if isconstrained(ctrl)
        updateQP_constrained!(ctrl, x, time)
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
    get_control(ctrl::MPCController, x, t)

Get the control from the MPC solver by solving the QP. 
If you want to use your own QP solver, you'll need to change this
method.
"""
function get_control(ctrl::MPCController{OSQP.Model}, x, time)
    # Update the QP
    updateQP!(ctrl, x, time)
    OSQP.update!(ctrl.solver, q=ctrl.q, l=ctrl.lb, u=ctrl.ub)

    # Solve QP
    results = OSQP.solve!(ctrl.solver)
    Δu = results.x[1:2]
    
    k = get_k(ctrl, time)
    return ctrl.Uref[k] + Δu 
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
    buildQP!(ctrl, A,B,Q,R,Qf; kwargs...)

Build the QP matrices `P` and `A` for the MPC problem. Note that these matrices
should be constant between MPC iterations.

Any keyword arguments will be passed to `initialize_solver!`.
"""
function buildQP_constrained!(ctrl::MPCController, A,B,Q,R,Qf; kwargs...)
    # TODO: Implement this method to build the QP matrices
    
    x_size = size(ctrl.Xref[1])[1]
    u_size = size(ctrl.Uref[1])[1]
    
    # useful variables
    N = ctrl.Nmpc
    R_size = size(R)
    Q_size = size(Q)
    P_rows = (N-1)*(R_size[1] + Q_size[1])
    P_cols = (N-1)*(R_size[2] + Q_size[2])
    
#     # populate first row of A matrix
#     A_constraint = [B -I(x_size) zeros(x_size, P_cols - u_size - x_size)]
    
#     # set first row of P
#     P = [R zeros(R_size[1], P_cols - R_size[2])]
#     curr_pos = u_size
    
#     # Inequality constraint matrices that will end up in A
#     theta_selector = zeros(x_size)
#     theta_selector[3] = 1
#     A_theta_constraint = [zeros(curr_pos)' theta_selector' zeros(P_cols - curr_pos - x_size)']
#     phi_selector = zeros(x_size)
#     phi_selector[8] = 1
#     A_phi_constraint = [zeros(curr_pos)' phi_selector' zeros(P_cols - curr_pos - x_size)']
#     z_selector = zeros(x_size)
#     z_selector[2] = 1
#     A_z_constraint = [zeros(curr_pos)' z_selector' zeros(P_cols - curr_pos - x_size)']
#     T_selector = zeros(x_size)
#     T_selector[7] = 1
#     A_T_constraint = [zeros(curr_pos)' T_selector' zeros(P_cols - curr_pos - x_size)']
    
#     for i=2:N-1
            
#         # Rows of A constraint matrix
#         new_A_const_row = [zeros(x_size, curr_pos) A B -I(x_size) zeros(x_size, P_cols - curr_pos - 2*x_size - u_size)]
#         A_constraint = [A_constraint; new_A_const_row]
        
#         # ==== P rows ====
#         # Q row
#         new_P_row = [zeros(Q_size[1], curr_pos) Q zeros(Q_size[1], P_cols - (curr_pos + Q_size[2]))]
#         curr_pos += x_size
#         P = [P; new_P_row]
            
#         # R row
#         new_P_row = [zeros(R_size[1], curr_pos) R zeros(R_size[1], P_cols - (curr_pos + R_size[2]))]
#         curr_pos += u_size
#         P = [P; new_P_row]
        
#         A_theta_constraint = [A_theta_constraint; zeros(curr_pos)' theta_selector' zeros(P_cols - curr_pos - x_size)']
#         A_phi_constraint = [A_phi_constraint; zeros(curr_pos)' phi_selector' zeros(P_cols - curr_pos - x_size)']
#         A_z_constraint = [A_z_constraint; zeros(curr_pos)' z_selector' zeros(P_cols - curr_pos - x_size)']
#         A_T_constraint = [A_T_constraint; zeros(curr_pos)' T_selector' zeros(P_cols - curr_pos - x_size)']
#     end
    
#     final_P_row = [zeros(Q_size[1], P_cols - Q_size[2]) Qf]
#     P = [P; final_P_row]
    
#     ctrl.P .= P
    
#     A_constraint = [A_constraint; A_theta_constraint; A_phi_constraint; A_z_constraint; A_T_constraint]
    
#     ctrl.A .= A_constraint
    
    
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
function updateQP_constrained!(ctrl::MPCController, x, time)
    
    # TODO: change this method to fit quadrotor
    
#     k = get_k(ctrl, time)
#     u_size = 2
#     x_size = 8
    
#     R = ctrl.P[begin:u_size, begin:u_size]
#     Q = ctrl.P[u_size+1:u_size+x_size, u_size+1:u_size+x_size]
#     Qf = ctrl.P[1+end-x_size:end, 1+end-x_size:end]
        
#     state_bounds = zeros(x_size*(ctrl.Nmpc-1))
#     state_bounds[begin:size(A)[1]] = -A*(x - xeq)
    
#     theta_ub = 5*pi/180.0 * ones(ctrl.Nmpc-1)
#     theta_lb = -theta_ub
#     phi_ub = 10*pi/180.0 * ones(ctrl.Nmpc-1)
#     phi_lb = -phi_ub
#     z_ub = Inf * ones(ctrl.Nmpc-1)
#     z_lb = zeros(ctrl.Nmpc-1)
#     T_ub = 0.5*model.m*model.g * ones(ctrl.Nmpc-1)
#     T_lb = -0.25*model.m*model.g * ones(ctrl.Nmpc-1)
    
#     ub = [state_bounds; theta_ub; phi_ub; z_ub; T_ub]
#     lb = [state_bounds; theta_lb; phi_lb; z_lb; T_lb]
                        
#     ctrl.ub .= ub
#     ctrl.lb .= lb
            
#     q = -R*(ctrl.Uref[k] - ueq)
            
#     # specify bounds for equality constraint
#     for i=1:ctrl.Nmpc-2
#         if i+k > size(ctrl.Xref)[1]
#             new_q_row = zeros(10)
#             q = [q; new_q_row]
#         else
#             new_q_row = -Q*(ctrl.Xref[k+i] - xeq)
#             q = [q; new_q_row]

#             new_q_row = -R*(ctrl.Uref[k+i] - ueq)
#             q = [q; new_q_row]
#         end
#     end
    
#     if k+ctrl.Nmpc-1 < size(ctrl.Xref)[1]
#         q = [q; -Qf*(ctrl.Xref[k+ctrl.Nmpc-1] - xeq)]
#     else
#         q = [q; zeros(8)]
#     end
                
#     ctrl.q .= q
                
    return nothing
end


"""
Create MPC object and solve
"""
function build_MPC_QP(Xref, Uref, tref, A, B, Q, R, Qf)
    # Initialize the constrained MPC controller
    Nd = (Nmpc-1)*(n+4)
    mpc = OSQPController(n, m, Nmpc, length(Xref), Nd)
    mpc.Xref .= Xref
    mpc.Uref .= Uref
    mpc.times .= tref
    
    # TODO: do we need to adjust A and B along the trajectory?
    buildQP!(mpc, A, B, Q, R, Qf, tol=1e-2, verbose=false)
    
    return mpc
end