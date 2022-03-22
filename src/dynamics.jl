# Taken from homework3-Q1

"""
Calculate the continuous time dynamics ẋ = f(x,u), x is a vector of length 6, u is a vector of length 2.

returns ẋ
"""

function dynamics(x,u,wind_disturbance=false)
    # planar quadrotor dynamics
    
    # parameters
    mass = 1.0 
    g = 9.81
    ℓ = 0.3 
    J = 0.2*mass*ℓ^2

     # unpack state
    px,pz,θ,vx,vz,ω = x
    
    # TODO: implement wind disturbance
    if wind_disturbance
        println("ERROR! wind_disturbance not implemented in dynamics")
        error("unimplemented")
        return zeros(6)
    else
        return [vx,vz,ω,(1/mass)*(u[1] + u[2])*sin(θ), (1/mass)*(u[1] + u[2])*cos(θ) - g, (ℓ/(2*J))*(u[2]-u[1])]
    end
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