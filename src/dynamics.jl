# Taken from homework3-Q1
using Random
using Setfield
using LinearAlgebra
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
returns ẋ
"""
function dynamics(x,u)
    # 3D quadrotor dynamics
    
    # unpack state
    r = x[1:3]
    q = x[4:7]
    v_B = x[8:10]
    ω_B = x[11:13]
    
    Q = rot_mat_from_quat(q)
    
    F_B = Q' * [0 0 -mass*g]' + [0 0 sum(k_T * u)]'    
    τ_B = τ_mat * u
    
    # see main.ipynb for formal dynamics definition
    ẋ = vcat(Q*v_B, 
             0.5*quat_L(q)*quat_H*ω_B,
             1/mass * F_B - cross(ω_B, v_B),
             inv(J)*(τ_B - cross(ω_B, J*ω_B)) )
    last_x = ẋ[:,1]
    return last_x
end

"""
Integrates the dynamics ODE 1 dt forward, x_{k+1} = rk4(x_k,u_k,dt).

returns x_{k+1}
"""
function rk4(x,u,dt, wind, wind_disturbance=false)
    # rk4 for integration
    k1 = dt*dynamics(x,u)
    k2 = dt*dynamics(x + k1/2,u)
    k3 = dt*dynamics(x + k2/2,u)
    k4 = dt*dynamics(x + k3,u)
    rk4_result = x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    if wind_disturbance
       rk4_result, win = simulate_random_walk_wind_traj!(rk4_result, wind)
       @set! wind.wd = win.wd
       @set! wind.wm = win.wm
       @set! wind.wind_dir = win.wind_dir
       return rk4_result, wind

    end
    return rk4_result
end

"""
uses forward diff to get the following jacobians of the above discrete dynamics function (rk4):
drk4/dx = A 
drk4/du = B
"""
function dynamics_jacobians(x,u,dt,wind, wind_disturbance=false )
    # returns the discrete time dynamics jacobians
    A = FD.jacobian(_x -> rk4(_x,u,dt, wind, wind_disturbance),x)
    B = FD.jacobian(_u -> rk4(x,_u,dt, wind, wind_disturbance),u)
    return A,B
end



function simulate_random_walk_wind_traj!(X, wind,move=0.1)
    rand_num = Random.rand()
    if rand_num < 0.5
       @set! wind.wd = wind.wd + move
    else
       @set! wind.wd = wind.wd - move
    end
    Fwind = RotMatrix{3}(RotZ(deg2rad(wind.wd))) * wind.wind_dir
    X[8] +=  Fwind[1]/mass 
    X[9] +=  Fwind[2]/mass 
    X[10] +=  Fwind[3]/mass 
    @set! wind.wind_dir = Fwind
    return X, wind
end    