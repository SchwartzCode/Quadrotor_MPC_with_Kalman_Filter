# Taken from homework3-Q1

"""
Quadrotor Dynamics Parameters
"""
mass = 1.0 # [kg]
g = 9.81   # gravitational acceleration [m/s^2]
ℓ = 0.3 # quadrotor arm length [m]
J = diagm([0.2, 0.2, 0.05]) * mass * ℓ^2 # quadrotor moments of inertia about (x,y,z)

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
function dynamics(x,u,wind_disturbance=false)
    # 3D quadrotor dynamics

     # unpack state
    r = x[1:3]
    q = x[4:7]
    v_B = x[8:10]
    ω_B = x[11:13]
    
    Q = rot_mat_from_quat(q)
    
#     println(Q)
    
#     println("arg", sum(k_T * u))
    F_B = Q' * [0 0 -mass*g]' + [0 0 sum(k_T * u)]'
    
#     println("\t", Q' * [0 0 -mass*g]')
#     println("\t- ", [0 0 sum(k_T * u)]')
    
#     println("= force: ", F_B)
    
    τ_B = τ_mat * u
    
    # see main.ipynb for formal dynamics definition
    ẋ = vcat(Q*v_B, 
             0.5*quat_L(q)*quat_H*ω_B,
             1/mass * F_B - cross(ω_B, v_B),
             inv(J)*(τ_B - cross(ω_B, J*ω_B)) )
    
    return ẋ
end

"""
Integrates the dynamics ODE 1 dt forward, x_{k+1} = rk4(x_k,u_k,dt).

returns x_{k+1}
"""
function rk4(x,u,dt)
    # rk4 for integration
    k1 = dt*dynamics(x,u)
    k2 = dt*dynamics(x + k1/2,u)
    k3 = dt*dynamics(x + k2/2,u)
    k4 = dt*dynamics(x + k3,u)
    return x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end

"""
uses forward diff to get the following jacobians of the above discrete dynamics function (rk4):
drk4/dx = A 
drk4/du = B
"""
function dynamics_jacobians(x,u,dt)
    # returns the discrete time dynamics jacobians
    A = FD.jacobian(_x -> rk4(_x,u,dt),x)
    B = FD.jacobian(_u -> rk4(x,_u,dt),u)
    return A,B
end