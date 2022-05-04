"""
    A collection of functions for help when using quaternions
    More info on these functions is available here: 
https://github.com/dynamics-simulation-16-715/lecture-notebooks/blob/main/Lecture%207/Lecture%207.pdf
    Which is the 7th lecture of the 2021 lectures for this course: https://github.com/dynamics-simulation-16-715/
"""

quat_H = zeros(4,3)
quat_H[2:4,:] = I(3)

function quat_L(q)
    [q[1] -q[2:4]'; q[2:4] q[1]*I + hat(q[2:4])]
end

function quat_R(q)
    [q[1] -q[2:4]'; q[2:4] q[1]*I - hat(q[2:4])]
end

function rot_mat_from_quat(q)
    # quat_H matrices used to trim off first row and column 
    #   (since R(q)' L(q) is 4x4 and Q is just bottom 3x3)
    return quat_H' * (quat_R(q)' * quat_L(q)) * quat_H
end

function get_quat_H()
    H = zeros(4,3)
    H[2:4,:] = I(3)
    return H
end

"""
Returns hat(ω)
"""
function hat(ω)
    return [0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end


function angle_of_rotation_q(q)
    angle = acos(q[1])*2
    return angle
end


function ρ(ϕ)
    # convert from ϕ to a quaternion 
    # read more: https://roboticexplorationlab.org/papers/planning_with_attitude.pdf
    return (1/sqrt(1 + dot(ϕ,ϕ)))*[1;ϕ]
end

function ϕ(q)
    # convert quaternion to rodrigues parameters
    
    return q[2:4] / q[1]
end

"""
Returns: attitude Jacobian given a quaternion
Params
    q: quaternion to calculate Jacobian about
"""
function G(q)
    return quat_L(q)*quat_H
end

"""
Returns: full attitude Jacobian matrix (i.e. padded with I to account
    for other state variables):
        Note: in this comment, 'x' refers to the length of state vector
         3  3  x-7
     3 | I  0   0 |
     4 | 0 G(Q) 0 |
   x-7 | 0  0   I |
Params
    x: 13-dim state vector (r, q, v_B, ω_B)
"""
function attitude_jacobian(x)
    x_size = length(x)
    
    q = x[4:7]
    
    return [I zeros(3,x_size-4); zeros(4,3) G(q) zeros(4,x_size-7); zeros(x_size-7,6) I]
end