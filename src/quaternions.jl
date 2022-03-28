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

function ρ(ϕ)
    # convert from ϕ to a quaternion 
    # read more: https://roboticexplorationlab.org/papers/planning_with_attitude.pdf
    return (1/sqrt(1 + dot(ϕ,ϕ)))*[1;ϕ]
end

"""
================================================================================================
Attitude Jacobian funcs below almost certainly need to be modified for our use:
================================================================================================
"""

function G(Q)
    return L(Q)*H
end

function Ḡ(q)
    Q = q[4:7]
    return [I zeros(3,3); zeros(4,3) G(Q)]
end