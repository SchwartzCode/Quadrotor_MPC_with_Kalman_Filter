# Taken from homework3-Q1
using Random
"""
Quadrotor Dynamics Parameters
"""
mass = 1.0 # [kg]
g = 9.81   # gravitational acceleration [m/s^2]
ℓ = 0.3 # quadrotor arm length [m]
J = diagm([0.004, 0.004, 0.008]) * mass * ℓ^2 # quadrotor moments of inertia about (x,y,z)

k_T = 1.0 #0.1 # [N/rpm]
k_m = 1.0 #0.1 # [N*m/rpm]
τ_mat = [[0 ℓ*k_T 0 -ℓ*k_T]
         [-ℓ*k_T 0 ℓ*k_T 0]
         [k_m -k_m k_m -k_m]]

"""
Calculate the continuous time dynamics ẋ = f(x,u), x is a vector of length 6, u is a vector of length 2.
param - x <13 elem vector>: [pos (world frame),
                             orientation (quaternion, body -> world),
                             velocity (body frame),
                             angular velocity (body frame)]
param - u <4 elem vector>: rpm of each motor (see main_script.ipynb for ordering)
param - wind_disturbance <boolean>: if true, wind disturbance is added to state derivative
returns ẋ

NOTE:
ANY CHANGES MADE TO THIS FUNCTION ALSO NEED TO BE MADE TO KF_dynamics_model()
IF THE KALMAN FILTER IS TO WORK
"""
function dynamics(x,u,wind)
    # 3D quadrotor dynamics
    
    # unpack state
    r = x[1:3]
    q = x[4:7]
    v_B = x[8:10]
    ω_B = x[11:13]
#     wind = x[14:16] # unused in this func
    
    Q = rot_mat_from_quat(q)
    F_B = Q' * ([0 0 -mass*g]' + reshape(wind.wind_dir, (3,1))) + [0 0 sum(k_T * u)]'
    τ_B = τ_mat * u
    
    # see main_script.ipynb for formal dynamics definition
    ẋ = vcat(Q*v_B, 
             0.5*quat_L(q)*quat_H*ω_B,
             1/mass * F_B - cross(ω_B, v_B),
             inv(J)*(τ_B - cross(ω_B, J*ω_B)),
             zeros(3)) # wind estimate is updated by Kalman filter; not here in dynamics    
    
    return ẋ[:,1]
end

"""
Same as dynamics() but uses state estimate of wind instead of actual wind
"""
function KF_dynamics_model(x,u)
    # 3D quadrotor dynamics
    
    # unpack state
    r = x[1:3]
    q = x[4:7]
    v_B = x[8:10]
    ω_B = x[11:13]
    wind = x[14:16]
    
    Q = rot_mat_from_quat(q)
    
    F_B = Q' * ([0 0 -mass*g]' + wind) + [0 0 sum(k_T * u)]'
    τ_B = τ_mat * u
    
    # see main_script.ipynb for formal dynamics definition
    ẋ = vcat(Q*v_B, 
             0.5*quat_L(q)*quat_H*ω_B,
             1/mass * F_B - cross(ω_B, v_B),
             inv(J)*(τ_B - cross(ω_B, J*ω_B)),
             zeros(3)) # FIXME: estimate of wind dynamics goes here
    
    return ẋ[:,1]
end

function simulate_random_walk_wind_traj!(wind, move=0.1)
    rand_num = Random.randn(1)[1]
    if wind.wind_type
        if rand_num > 0.5
            wind.wind_dir .= RotMatrix{3}(RotZ(move)) * wind.wind_dir
        else
            wind.wind_dir .= RotMatrix{3}(RotZ(-move)) * wind.wind_dir
        end
    else
        if rand_num < 0.5
           wind.wd[1] = wind.wd[1] + move
        else
           wind.wd[1] = wind.wd[1] - move
        end
        wind.wind_dir .= RotMatrix{3}(RotZ(deg2rad(wind.wd[1]))) * wind.wind_dir
    end
end   

function get_wind_correction(x,B)
    # unpack state
    r = x[1:3]
    q = x[4:7]
    v_B = x[8:10]
    ω_B = x[11:13]
    wind = x[14:16]
    
    Q = rot_mat_from_quat(q)
    F_wind = Q' * wind # body frame
    
    dv_du = B[8:10,:] # 3x4
    
    thrust_correction = dv_du \ -F_wind
    
#     println("wind force: ", wind)
#     println("thrust corr: ", thrust_correction)
    
    return thrust_correction
end

"""
Integrates the dynamics ODE 1 dt forward, x_{k+1} = rk4(x_k,u_k,dt).

returns x_{k+1}
"""
function rk4(x,u,dt,wind)
    # rk4 for integration
    k1 = dt*dynamics(x,u, wind)
    k2 = dt*dynamics(x + k1/2,u, wind)
    k3 = dt*dynamics(x + k2/2,u, wind)
    k4 = dt*dynamics(x + k3,u, wind)
    return x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end

"""
uses forward diff to get the following jacobians of the above discrete dynamics function (rk4):
drk4/dx = A 
drk4/du = B
"""
function dynamics_jacobians(x,u,dt,wind)
    # returns the discrete time dynamics jacobians
    A = FD.jacobian(_x -> rk4(_x,u,dt,wind),x)
    B = FD.jacobian(_u -> rk4(x,_u,dt,wind),u)
    return A,B
end


"""
Integrates the dynamics ODE 1 dt forward, x_{k+1} = rk4(x_k,u_k,dt).

returns x_{k+1}
"""
function KF_rk4(x,u,dt)
    # rk4 for integration
    k1 = dt*KF_dynamics_model(x,u)
    k2 = dt*KF_dynamics_model(x + k1/2,u)
    k3 = dt*KF_dynamics_model(x + k2/2,u)
    k4 = dt*KF_dynamics_model(x + k3,u)
    return x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end


function KF_dynamics_jacobians(x,u,dt)
    # returns the discrete time dynamics jacobians using state wind estimate
    A = FD.jacobian(_x -> KF_rk4(_x,u,dt),x)
    B = FD.jacobian(_u -> KF_rk4(x,_u,dt),u)
    return A,B
end
 