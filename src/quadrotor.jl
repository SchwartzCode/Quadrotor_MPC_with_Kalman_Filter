using MeshCat
using RobotZoo: Quadrotor, PlanarQuadrotor
using CoordinateTransformations, Rotations, Colors, StaticArrays, RobotDynamics

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
    n = state_dim(model)
    for (k,x) in enumerate(X)
        atframe(anim, k) do
            x = X[k]
            visualize!(vis, model, SVector{n}(x)) 
        end
    end
    setanimation!(vis, anim)
end

# TODO: re-work this func to be 3D (not planar)
"""
params
N  - number of time steps
"""
function line_reference(N::Int64, dt::Float64)
    x_ref = Array(zeros(N))
    y_ref = Array(LinRange(-3,3,N))
    z_ref = Array(ones(N))
    quat_ref = hcat(ones(N), zeros(N), zeros(N), ones(N))

    vx_ref = Array(zeros(N))
    vy_ref = Array((6.0/(N*dt))*ones(N))
    vz_ref = Array(zeros(N))
    ωx_ref = Array(zeros(N))
    ωy_ref = Array(zeros(N))
    ωz_ref = Array(zeros(N))
    
    xref = [x_ref'; y_ref'; z_ref'; quat_ref'; vx_ref'; vy_ref'; vz_ref'; ωx_ref'; ωy_ref'; ωz_ref'] 
    return [SVector{13}(x) for x in eachcol(xref)]
end

"""
params
N  - number of time steps
dt - size of time steps
"""
# TODO: flip looks wonky
function flip_reference(N::Int64, dt::Float64)
    N_pre_flip = Int(floor(N / 4))
    N_flip = Int(floor(N / 2))
    N_post_flip = Int(ceil(N / 4))
    N_post_flip += N - (N_pre_flip + N_flip + N_post_flip)
    half_N_flip = Int(N_flip/2)
    
    x_ref = Array(zeros(N))
    y_ref = [LinRange(-3,0,N_pre_flip); zeros(N_flip); LinRange(0,3,N_post_flip)]
    z_ref = [ones(N_pre_flip); LinRange(1,3,half_N_flip); LinRange(3,1,half_N_flip); ones(N_post_flip)]
    
    # TODO: populate this properly
    quat_ref = hcat(ones(N), zeros(N), zeros(N), ones(N))
    println(N_flip)
    flip_angles = collect(LinRange(0.0, -2*pi, N_flip))
    
    quat_ref[N_pre_flip+1:N_pre_flip+N_flip,1] = cos.(flip_angles/2)
    quat_ref[N_pre_flip+1:N_pre_flip+N_flip,3] = sin.(flip_angles/2)
    
    vx_ref = Array(zeros(N))
    vy_ref = [3.0/(N_pre_flip/dt)*ones(N_pre_flip); zeros(N_flip); 3.0*ones(N_post_flip)]
    vz_ref = [zeros(N_pre_flip); 2.0/(half_N_flip/dt)*ones(half_N_flip); -2.0/half_N_flip/dt*ones(half_N_flip); zeros(N_post_flip)]
    ωy_ref = [zeros(N_pre_flip); -2*pi/(N_flip/dt)*ones(N_flip); zeros(N_post_flip)]
    ωx_ref = Array(zeros(N))
    ωz_ref = Array(zeros(N))
    
    xref = [x_ref'; y_ref'; z_ref'; quat_ref'; vx_ref'; vy_ref'; vz_ref'; ωx_ref'; ωy_ref'; ωz_ref'] 
    return [SVector{13}(x) for x in eachcol(xref)]
end

function RobotDynamics.discrete_jacobian!(::Type{Q}, ∇f, model::AbstractModel,
        x, u, t, dt) where {Q<:RobotDynamics.Explicit}
    z = KnotPoint(x, u, dt, t)
    RobotDynamics.discrete_jacobian!(Q, ∇f, model, z)
end

struct WindyQuad <: AbstractModel
    quad::Quadrotor
    dir::MVector{2,Float64}   # wind direction
    wd::Float64               # std on wind angle
    wm::Float64               # std on wind magnitude
end

# TODO @Corinne - will need to change this struct to be 3D
function WindyQuad(quad::Quadrotor;
        wind = [1,1]*1.0,
        wd = deg2rad(10),
        wm = 0.01,
    )
    WindyQuad(quad, SA_F64[wind[1], wind[2]], Float64(wd), Float64(wm)) 
end
RobotDynamics.state_dim(model::WindyQuad) = state_dim(model.quad)
RobotDynamics.control_dim(model::WindyQuad) = control_dim(model.quad)

# TODO @Corinne - will need to re-work this function to be 3D (not planar)
function RobotDynamics.dynamics(model::WindyQuad, x, u)
    ẋ = dynamics(model.quad, x, u)
    mass = model.quad.mass
    wind_mag = randn()*model.wm
    wind_dir = Angle2d(randn()*model.wd) * model.dir
    Fwind =  Angle2d(randn()*model.wd) * model.dir 
    ẋ2 = SA[ẋ[1], ẋ[2], ẋ[3], ẋ[4] + Fwind[1]/mass, ẋ[5] + Fwind[2]/mass, ẋ[6]]
    return ẋ2
end

function simulate(quad::Quadrotor, x0, ctrl; tf=1.5, dt=0.025, kwargs...)
    model = WindyQuad(quad; kwargs...)

    n,m = size(model)
    times = range(0, tf, step=dt)
    N = length(times)
    X = [@SVector zeros(n) for k = 1:N] 
    U = [@SVector zeros(m) for k = 1:N-1]
    X[1] = x0

    tstart = time_ns()

    for k = 1:N-1
        U[k] = get_control(ctrl, X[k], times[k])
        X[k+1] = discrete_dynamics(RK4, model, X[k], U[k], times[k], dt)
    end
    tend = time_ns()
    rate = N / (tend - tstart) * 1e9
    println("Controller ran at $rate Hz")
    return X,U,times
end
