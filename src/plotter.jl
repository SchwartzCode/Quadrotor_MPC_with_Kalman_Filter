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
        xlabel="Time [sec]", ylabel="Position [m]")
    plot!(tref, vals_act, label=act_labels)
end