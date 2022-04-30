function plot_vals(X_ref, X_act, indices, labels, plot_title="title")
    
    dims = length(collect(indices))
    # extract values to plot from given matrices
    vals_ref = zeros(N,dims)
    vals_ref[1,:] .= X_ref[1][indices]
    vals_act = zeros(N,dims)
    vals_act[1,:] .= X_act[1][indices]

    for k = 1:N-1
        vals_ref[k+1,:] .= X_ref[k+1][indices]
        vals_act[k+1,:] .= X_act[k+1][indices]
    end
    
    # initialize label mats and populate (there's probably a better way to do this)
    ref_labels = Array{Union{Nothing, String}}(nothing, 1, length(labels))
    act_labels = Array{Union{Nothing, String}}(nothing, 1, length(labels))
    for i=1:length(labels)
        ref_labels[i] = labels[i]*"_ref"
        act_labels[i] = labels[i]*"_act"
    end
    
    plot(tref, vals_ref, line=:dash, label=ref_labels, title=plot_title,
        xlabel="Time [sec]", ylabel="Value")
    plot!(tref, vals_act, label=act_labels)
end


function plot_wind_tracking(X_KF, wind_vals, plot_title="Wind Force Estimation")
    N = length(wind_hist[:,1])
    
    # extract values to plot from given matrices
    vals_est = zeros(N,3)
    vals_est[1,:] .= X_KF[1][14:16]
    
    for k = 1:N-1
        vals_est[k+1,:] .= X_KF[k+1][14:16]
    end
    
    # initialize label mats and populate (there's probably a better way to do this)
    est_labels = ["wind_x est" "wind_y est" "wind_z est"]
    gt_labels = ["wind_x actual" "wind_y actual" "wind_z actual"]
    
    plot(tref, wind_hist, line=:dash, label=gt_labels, title=plot_title,
        xlabel="Time [sec]", ylabel="Force [N]", legend=:bottomleft)
    plot!(tref, vals_est, label=est_labels)
end