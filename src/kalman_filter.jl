"""
I used slide 63 of Prof Kaess' EKF notes as a guide when writing this EKF code
https://canvas.cmu.edu/courses/27971/files/7701028?module_item_id=5089860
"""

"""
x/u/Σ are from current state, want to predict x and Σ at next state
"""
function EKF_predict(x, u, Σ, dt)
    x_next = KF_rk4(x,u,dt)
    A,B = KF_dynamics_jacobians(x,u,dt)
    J_att = attitude_jacobian(x)
    A = J_att'*A*J_att
    
    Σ_next = A*Σ*A'

    return x_next, Σ_next
end

# note: will need to add measurement noise to readings, otherwise covariance matrix becomes singular
Q_t = 1e-12*I(15)
# Q_t[7:9,7:9] = 0.1*I(3)
# Q_t[13:15,13:15] = 5*I(3)

"""
correct estimate given error in state vs actual state
"""
function EKF_correct(x_pred, u, measurement, Σ_pred, dt)
    
    # measurement jacobian
    H_t = [I(12) zeros(12,3); zeros(3,15)]
    
    # update measurement jacobian with effect wind has on velocity
    A,B = KF_dynamics_jacobians(x_pred,u,dt)
    J_att = attitude_jacobian(x_pred)
    A = J_att'*A*J_att
    H_t[7:9,13:15] .= -A[7:9,13:15]
    println(A[7:9,13:15])
#     H_t = [zeros(3,7) A[8:10,14:16] zeros(3,6)]
    
#     println(size(H_t))
#     println(size(Σ_pred))
    K = Σ_pred*H_t'*inv(H_t*Σ_pred*H_t' + Q_t)
    
    err = EKF_err(measurement, x_pred)
#     println("KF update: ", K * err)
#     println("diff: ", err)
    
    # measurement is true state except for wind
#     update = vcat([(K * err)[1:3]' ρ((K*err)[4:6])' (K * err)[7:end]'])'
#     println("\n update: ", update)
#     x_corrected = x_pred + update[:,1]
#     wind_update = x_pred + K * (measurement[8:10] - x_pred[8:10])
    diff = (measurement - x_pred)
#     diff = setindex(diff, 4, 0.0)
#     diff[4:7] = zeros(4)#0.0
        wind_update = x_pred[2:end] + K * err #FIXME: dirty
    println(norm(K[14:end,:]))
    x_corrected = vcat([measurement[1:13]' wind_update[13:15]'])'[:,1]
    
    Σ_corrected = (I - K*H_t)*Σ_pred + K*Q_t*K'
    
    return x_corrected, Σ_corrected
end

"""
difference between state but converts quaternions to 3 param
attitude representation using Cayley Map
"""
function EKF_err(measurement, x_pred)
    dϕ = ϕ(quat_L(x_pred[4:7])' * measurement[4:7])
    
    dX = vcat([(measurement - x_pred)[1:3]' dϕ' (measurement - x_pred)[8:end]'])'
    
    return dX
end