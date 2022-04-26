using MeshCat
using RobotZoo: Quadrotor, PlanarQuadrotor
using CoordinateTransformations, Rotations, Colors, StaticArrays, RobotDynamics
using Random
using TrajectoryOptimization
# using TrajOptPlots
using Altro

function set_mesh!(vis, model::L;
        scaling=1.0, color=colorant"black"
    ) where {L <: Union{Quadrotor, PlanarQuadrotor}} 
    # urdf_folder = joinpath(@__DIR__, "..", "data", "meshes")
    urdf_folder = @__DIR__
    # if scaling != 1.0
    #     quad_scaling = 0.085 * scaling
    obj = joinpath(urdf_folder, "quadrotor_scaled.obj")
    if scaling != 1.0
        error("Scaling not implemented after switching to MeshCat 0.12")
    end
    robot_obj = MeshFileGeometry(obj)
    mat = MeshPhongMaterial(color=color)
    setobject!(vis["robot"]["geom"], robot_obj, mat)
    if hasfield(L, :ned)
        model.ned && settransform!(vis["robot"]["geom"], LinearMap(RotX(pi)))
    end
end

function visualize!(vis, model::Quadrotor, x::StaticVector)
    px,py,pz = x[1], x[2], x[3]
    rot_mat = rot_mat_from_quat(x[4:7])
    
    settransform!(vis["robot"], compose(Translation(px,py,pz), LinearMap(rot_mat)))
end

function visualize!(vis, model, tf::Real, X)
    fps = Int(round((length(X)-1)/tf))
    anim = MeshCat.Animation(fps)
    n,m = RobotDynamics.dims(model)
    for (k,x) in enumerate(X)
        atframe(anim, k) do
            x = X[k]
            visualize!(vis, model, SVector{n}(x)) 
        end
    end
    setanimation!(vis, anim)
end



function altro_reference_line(N::Int64, dt::Float64, model)
    n,m = RobotDynamics.dims(model)
    # Define initial and final conditions
    x0 = @SVector [0, -1, 0.5, 1, 0, 0, 0, 2.0/(N*dt), 2.0/(N*dt),2.0/(N*dt),0, 0, 0] 
    xf = @SVector [2, 1, 2.5, 1, 0, 0, 0, 2.0/(N*dt), 2.0/(N*dt),2.0/(N*dt),0, 0, 0] 
    Q = 1.0e-2*Diagonal(@SVector ones(n)) * dt
    Qf = 100.0*Diagonal(@SVector ones(n))
    R = 1.0e-1*Diagonal(@SVector ones(m)) * dt
    println(size(Q))
    println(size(xf))
    # Set up
    obj = LQRObjective(Q,R,Qf,xf,N)
    # Add constraints
    conSet = ConstraintList(n,m,N)
    u_bnd = 3.0
        # TODO: move this to build
    thrust_ub = 10.0
    thrust_lb = -5.0
  
    bnd = BoundConstraint(n, m, u_min=thrust_lb, u_max=thrust_ub)
    goal = GoalConstraint(xf)
    add_constraint!(conSet, bnd, 1:N-1)
    add_constraint!(conSet, goal, N)
    # Initialization
    u0 = @SVector fill(0.01,m)
    U0 = [u0 for k = 1:N-1]
    # Define problem
    tf = (N-1)* dt
    prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
    initial_controls!(prob, U0)
    # Solve with ALTRO
    opts = SolverOptions(
        cost_tolerance_intermediate=1e-2,
        penalty_scaling=10.,
        penalty_initial=1.0
    )
    altro = ALTROSolver(prob, opts)
    solve!(altro)
    # Extract the solution
    X = states(altro)
    println(X)
    return X
end
function altro_reference_circle(N::Int64, dt::Float64, model)
    N_pre_flip = Int(floor(N / 4))
    N_flip = Int(floor(N / 2 /4))
    println(N_flip)
    N_post_flip = Int(ceil(N / 4))
    N_post_flip += N - (N_pre_flip + N_flip + N_post_flip)
   
        
#     vx_ref = Array(zeros(N))
#     vy_ref = [vy_1; zeros(N_flip); vy_2]
#     vz_ref = [zeros(N_pre_flip); vz_1; vz_2; zeros(N_post_flip)]
#     ωx_ref = [zeros(N_pre_flip); ωx_flip; zeros(N_post_flip)]
#     ωy_ref = Array(zeros(N))
#     ωz_ref = Array(zeros(N))
    # Define initial and final conditions
    #fly to position to start flip
                 #x  y    z  w  x  y  z  vx  vy        vz ωx ωy ωz        
    x0 = @SVector [0, -3, 1, 1, 0, 0, 0, 0, 2.0/(N*dt), 0, 0, 0, 0] 
    xf = @SVector [0, 0,  1, 1, 0, 0, 0, 0, 2.0/(N*dt), 0, 0, 0, 0] 
    X1 = altro_reference(N_pre_flip, dt, model, x0, xf)
    #fly to pos 1 of flip
                 #x  y    z  w  x  y  z  vx  vy        vz ωx ωy ωz        
    x0 = X1[end]
    xf = @SVector [0, 3,  3, 1, 0, 0, 0, 0, 2.0/(N*dt),2.0/(N*dt), 0, 0, 0] 
    X2 = altro_reference(N_flip, dt, model, x0, xf)
    
#     #fly to top of flip    
                 #x  y    z  w  x  y  z  vx  vy        vz ωx ωy ωz        
    x0 = X2[end]
    xf = @SVector [0, 0,  6, 1, 0, 0, 0, 0, 2.0/(N*dt),2.0/(N*dt), 0, 0, 0] 
    X3 = altro_reference(N_flip, dt, model, x0, xf)
    
   #fly to pos2 of flip    
                  #x  y    z  w  x  y  z  vx  vy        vz ωx ωy ωz        
    x0 = X3[end]
    xf = @SVector [0, -3,  3, 1, 0, 0, 0, 0, 2.0/(N*dt),2.0/(N*dt), 0, 0, 0] 
    X4 = altro_reference(N_flip, dt, model, x0, xf)
    
    #fly to end of flip    
                  #x  y    z  w  x  y  z  vx  vy        vz ωx ωy ωz        
    x0 = X4[end]
    xf = @SVector [0, 0,  1, 1, 0, 0, 0, 0, 2.0/(N*dt),2.0/(N*dt), 0, 0, 0] 
    X5 = altro_reference(N_flip, dt, model, x0, xf)
    
    #fly to goal
                     #x  y    z  w  x  y  z  vx  vy        vz ωx ωy ωz        
    x0 = X5[end] 
    xf = @SVector [0, 3,  1, 1, 0, 0, 0, 0, 2.0/(N*dt),2.0/(N*dt), 0, 0, 0] 
    X6 = altro_reference(N_post_flip, dt, model, x0, xf)
    
    X = vcat(X1, X2,X3, X4,X5, X6)
    println(size(X))
    return X
    
#                   #x  y   z    w  x  y  z  xv          vy          vz          ωx ωy ωz        
#     x0 = @SVector [0, -1, 0.5, 1, 0, 0, 0, 2.0/(N*dt), 2.0/(N*dt), 2.0/(N*dt), 0, 0, 0] 
#     xf = @SVector [2, 1,  2.5, 1, 0, 0, 0, 2.0/(N*dt), 2.0/(N*dt), 2.0/(N*dt), 0, 0, 0] 
end
    
 
function altro_reference(N::Int64, dt::Float64, model, x0, xf)
    n,m = RobotDynamics.dims(model)
    Q = 1.0e-2*Diagonal(@SVector ones(n)) * dt
    Qf = 100.0*Diagonal(@SVector ones(n))
    R = 1.0e-1*Diagonal(@SVector ones(m)) * dt
    println(size(Q))
    println(size(xf))
    # Set up
    obj = LQRObjective(Q,R,Qf,xf,N)
    # Add constraints
    conSet = ConstraintList(n,m,N)
        
    #TODO: MAKE THIS NOT MAGIC NUMBERS
    u_bnd = 3.0
    thrust_ub = 10.0
    thrust_lb = -5.0
  
    bnd = BoundConstraint(n, m, u_min=thrust_lb, u_max=thrust_ub)
    goal = GoalConstraint(xf)
    add_constraint!(conSet, bnd, 1:N-1)
    add_constraint!(conSet, goal, N)
    # Initialization
    u0 = @SVector fill(0.01,m)
    U0 = [u0 for k = 1:N-1]
    # Define problem
    tf = (N-1)* dt
    prob = Problem(model, obj, x0, tf, xf=xf, constraints=conSet)
    initial_controls!(prob, U0)
    # Solve with ALTRO
    opts = SolverOptions(
        cost_tolerance_intermediate=1e-2,
        penalty_scaling=10.,
        penalty_initial=1.0
    )
    altro = ALTROSolver(prob, opts)
    solve!(altro)
    # Extract the solution
    X = states(altro)
    return X
end






function trapezoidal_vel(N::Int64, dt::Float64, p0, pf)
    quarter_N = Int(floor(N / 4))
    half_N = N - 2*quarter_N
    
    V_avg = 4*(pf - p0)/(3*N*dt)
    
    vels = [LinRange(0,V_avg,quarter_N); V_avg*ones(half_N); LinRange(V_avg,0,quarter_N)]
    positions = cumsum(vels, dims=1)*dt .+ p0
    
    return positions, vels
end

"""
params
N  - number of time steps
"""
function line_reference(N::Int64, dt::Float64)
    x_ref = Array(LinRange(0,2,N))
    y_ref = Array(LinRange(-1,1,N))
    z_ref = Array(LinRange(0.5,2.5,N))
    quat_ref = hcat(ones(N), zeros(N), zeros(N), zeros(N))

    vx_ref = Array((2.0/(N*dt))*ones(N))
    vy_ref = Array((2.0/(N*dt))*ones(N))
    vz_ref = Array((2.0/(N*dt))*ones(N))
    ωx_ref = Array(zeros(N))
    ωy_ref = Array(zeros(N))
    ωz_ref = Array(zeros(N))
    
    xref = [x_ref'; y_ref'; z_ref'; quat_ref'; vx_ref'; vy_ref'; vz_ref'; ωx_ref'; ωy_ref'; ωz_ref'] 
    return [SVector{13}(x) for x in eachcol(xref)]
end

function circle_trajectory(z,r)
    N = length(z)
    y = zeros(N)
    half_N = trunc(Int, N/4)
    y_1 = zeros(half_N)
   
    for i=1:half_N
        y_1[i] = real((Complex(r^2 - (z[i])^2))^0.5)
    end
    y[1:half_N] = reverse(y_1)
    y[half_N+1:2*half_N] = y_1
    y[2*half_N+1:3*half_N] = -reverse(y_1)
    y[3*half_N+1:N] = -y_1
    return y
    
end
"""
params
N  - number of time steps
dt - size of time steps
"""
function flip_reference(N::Int64, dt::Float64)
    # TODO: make pre and post flip N the same and add a tail to traj
    N_pre_flip = Int(floor(N / 4))
    N_flip = Int(floor(N / 2))
    N_post_flip = Int(ceil(N / 4))
    N_post_flip += N - (N_pre_flip + N_flip + N_post_flip)
    half_N_flip = Int(N_flip/2)
    
    # trapezoidal vels for y pos/vel components
    py_1, vy_1 = trapezoidal_vel(N_pre_flip, dt, -3, 0)
#     py_1 = py_1 .- 3
    py_2, vy_2 = trapezoidal_vel(N_post_flip, dt, 0, 3)
    
    # trapezoidal vels for z pos/vel components
    pz_1, vz_1 = trapezoidal_vel(half_N_flip, dt, 1, 3)
#     pz_2, vz_2 = trapezoidal_vel(half_N_flip, dt, 3, 1)
#     pz_1 .+ 1.0
    pz_2 = reverse(pz_1)
    vz_2 = -vz_1
    
    x_ref = Array(zeros(N))
    z_flip = [LinRange(1,3,half_N_flip); LinRange(3,1,half_N_flip)]
    r = 2
    y_flip = circle_trajectory(z_flip,r);
    y_ref = [LinRange(-3,y_flip[1],N_pre_flip); y_flip; LinRange(y_flip[end],3,N_post_flip)]
    z_ref = [ones(N_pre_flip); z_flip; ones(N_post_flip)]
    
    # init with quaternion of all 0 rad euler angles
    quat_ref = hcat(ones(N), zeros(N), zeros(N), zeros(N))
    
    # do full rotation about y axis
#     angs = collect(LinRange(0, 2*pi, N_flip))
    angs, ωx_flip = trapezoidal_vel(N_flip, dt, 0, 2*pi)
    flip_angles = [tan.(angs/2) zeros(N_flip) zeros(N_flip)]
    
    for i=1:N_flip
        quat_ref[N_pre_flip+i,:] .= ρ(flip_angles[i,:])
    end
    
    vx_ref = Array(zeros(N))
    vy_ref = [vy_1; zeros(N_flip); vy_2]
    vz_ref = [zeros(N_pre_flip); vz_1; vz_2; zeros(N_post_flip)]
    ωx_ref = [zeros(N_pre_flip); ωx_flip; zeros(N_post_flip)]
    ωy_ref = Array(zeros(N))
    ωz_ref = Array(zeros(N))
    
    vels = [vx_ref'; vy_ref'; vz_ref']
    
    for i = 1:N
        Q = rot_mat_from_quat(quat_ref[i,:])
        vels[:,i] = Q * vels[:,i]
    end
    
    xref = [x_ref'; y_ref'; z_ref'; quat_ref'; vels; ωx_ref'; ωy_ref'; ωz_ref'] 
    return [SVector{13}(x) for x in eachcol(xref)]
end



"""
params
N  - number of time steps
dt - size of time steps
"""
function other_flip_reference(N::Int64, dt::Float64)
    # TODO: make pre and post flip N the same and add a tail to traj
    N_pre_flip = Int(floor(N / 4))
    N_flip = Int(floor(N / 2))
    N_post_flip = Int(ceil(N / 4))
    N_post_flip += N - (N_pre_flip + N_flip + N_post_flip)
    half_N_flip = Int(N_flip/2)
    
    # trapezoidal vels for y pos/vel components
    py_1, vy_1 = trapezoidal_vel(N_pre_flip, dt, -3, 0)
#     py_1 = py_1 .- 3
    py_2, vy_2 = trapezoidal_vel(N_post_flip, dt, 0, 3)
    
    # trapezoidal vels for z pos/vel components
    pz_1, vz_1 = trapezoidal_vel(half_N_flip, dt, 1, 3)
#     pz_2, vz_2 = trapezoidal_vel(half_N_flip, dt, 3, 1)
#     pz_1 .+ 1.0
    pz_2 = reverse(pz_1)
    vz_2 = -vz_1
    
    x_ref = Array(zeros(N))
    y_ref = [py_1; zeros(N_flip); py_2]
    z_ref = [ones(N_pre_flip); pz_1; pz_2; ones(N_post_flip)]
    
    # init with quaternion of all 0 rad euler angles
    quat_ref = hcat(ones(N), zeros(N), zeros(N), zeros(N))
    
    # do full rotation about y axis
#     angs = collect(LinRange(0, 2*pi, N_flip))
    angs, ωx_flip = trapezoidal_vel(N_flip, dt, 0, 2*pi)
    flip_angles = [tan.(angs/2) zeros(N_flip) zeros(N_flip)]
    
    for i=1:N_flip
        quat_ref[N_pre_flip+i,:] .= ρ(flip_angles[i,:])
    end
    
    vx_ref = Array(zeros(N))
    vy_ref = [vy_1; zeros(N_flip); vy_2]
    vz_ref = [zeros(N_pre_flip); vz_1; vz_2; zeros(N_post_flip)]
    ωx_ref = [zeros(N_pre_flip); ωx_flip; zeros(N_post_flip)]
    ωy_ref = Array(zeros(N))
    ωz_ref = Array(zeros(N))
    
    vels = [vx_ref'; vy_ref'; vz_ref']
    
    for i = 1:N
        Q = rot_mat_from_quat(quat_ref[i,:])
        vels[:,i] = Q * vels[:,i]
    end
    
    xref = [x_ref'; y_ref'; z_ref'; quat_ref'; vels; ωx_ref'; ωy_ref'; ωz_ref'] 
    return [SVector{13}(x) for x in eachcol(xref)]
end

# function RobotDynamics.discrete_jacobian!(::Type{Q}, ∇f, model::AbstractModel,
#         x, u, t, dt) where {Q<:RobotDynamics.Explicit}
#     z = KnotPoint(x, u, dt, t)
#     RobotDynamics.discrete_jacobian!(Q, ∇f, model, z)
# end

# struct WindyQuad <: AbstractModel
#     quad::Quadrotor
#     dir::MVector{2,Float64}   # wind direction
#     wd::Float64               # std on wind angle
#     wm::Float64               # std on wind magnitude
# end

# # TODO @Corinne - will need to change this struct to be 3D
# function WindyQuad(quad::Quadrotor;
#         wind = [1,1]*1.0,
#         wd = deg2rad(10),
#         wm = 0.01,
#     )
#     WindyQuad(quad, SA_F64[wind[1], wind[2]], Float64(wd), Float64(wm)) 
# end
# RobotDynamics.state_dim(model::WindyQuad) = state_dim(model.quad)
# RobotDynamics.control_dim(model::WindyQuad) = control_dim(model.quad)

function simulate(quad::Quadrotor, x0, ctrl, A, B; tf=1.5, dt=0.025, kwargs...)
#     model = WindyQuad(quad; kwargs...)

    n,m = size(quad)
    times = range(0, tf, step=dt)
    N = length(times)
    X = [@SVector zeros(n) for k = 1:N] 
    U = [@SVector zeros(m) for k = 1:N-1]
    X[1] = x0

    tstart = time_ns()

    for k = 1:N-1
        U[k] = get_control(ctrl, A, B, X[k], times[k])
        X[k+1] = rk4(X[k], U[k], dt)
    end
    tend = time_ns()
    rate = N / (tend - tstart) * 1e9
    println("Controller ran at $rate Hz")
    return X,U,times
end

