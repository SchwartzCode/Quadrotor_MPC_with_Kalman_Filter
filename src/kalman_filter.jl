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
    Σ_next = A*Σ*A'

    return x_next, Σ_next
end

# note: will need to add measurement noise to readings, otherwise covariance matrix becomes singular
Q_t = 1e-6*I(16)
# Q_t[14:16,14:16] = 0.1*I(3)

"""
correct estimate given error in state vs actual state
"""
function EKF_correct(x_pred, measurement, Σ_pred, dt)
    
    H_t = [I(13) zeros(13,3); zeros(3,16)]
    # effect wind has on velocity
    Q = rot_mat_from_quat(measurement[4:7])
    H_t[8:10,14:16] = Q'
    
    K = Σ_pred*H_t'*inv(H_t*Σ_pred*H_t' + Q_t)
    
    println("Update: ", K * (measurement - x_pred))
    # measurement is true state, skip 
    x_corrected = x_pred + K * (measurement - x_pred)
    
    Σ_corrected = (I - K*H_t)*Σ_pred
    
    return x_corrected, Σ_corrected
end