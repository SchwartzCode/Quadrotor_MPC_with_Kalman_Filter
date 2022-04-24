
using MeshCat
using CoordinateTransformations
using Rotations 
using Colors
using GeometryBasics
using RobotZoo: PlanarRocket
using Plots


#taken from hw2.q3
function comparison_plot(model, Z...)
    #usage:  comparison_plot(model, (Xlqr,Ulqr,tlqr,"LQR"), (Xmpc1,Umpc1,tmpc1,"MPC"))
    p = plot(layout=(2,2), size=(1000,800))
    for z in Z
        plot!(p[1], z[3],  z[1],  inds=1:3, label=["x", "y", "z"],
            xlabel="time (s)", ylabel="Positions", legend=:topright)
        plot!(p[2], z[3],  z[1], inds=4:7, label=["q0", "q1", "q2", "q3"],
            xlabel="time (s)", ylabel="attitude(deg)", legend=:none)
        plot!(p[3], z[3],  z[1], inds=8:10, label=z[4],
            xlabel="time (s)", ylabel="Thrust (N)", legend=:none)
        plot!(p[4], z[3],  z[1], inds=11:13, label=["ωx", "ωy", "ωz"],
            xlabel="time (s)", ylabel="Thrust angle (deg)", legend=:none)
    end
   p 
end

function plot_trajectories(ref_traj, act_traj)
    plt = plot3d(
        2,
        xlim = (-3, 3),
        ylim = (-3, 3),
        zlim = (-5, 6),
        title = "Trajectory",
        marker = 2,
    )
#     plot!(plt, ref_traj[:][1],ref_traj[:][2], ref_traj[:][3],label="test")
#     plot!(plt, act_traj[:][1],act_traj[:][2], act_traj[:][3])
    
    # build an animated gif by pushing new points to the plot, saving every 10th frame
    @gif for i=1:N
#         println(ref_traj[i][1]," ",ref_traj[i][2], " ",ref_traj[i][3])
        plot!(plt, ref_traj[i][1],ref_traj[i][2], ref_traj[i][3],label="test")
        plot!(plt, act_traj[i][1],act_traj[i][2], act_traj[i][3])
    end every 10
#     plot(tref, vals_ref, line=:dash, label=["x_ref" "y_ref" "z_ref"], title="Position Plots",
#     xlabel="Time [sec]", ylabel="Position [m]")
#     plot!(tref, vals_act, label=["x_act" "y_act" "z_act"])
    
    
#     plot(tref, vals_ref, line=:dash, label=["x_ref" "y_ref" "z_ref"], title="Velocity Plots",
#     xlabel="Time [sec]", ylabel="Position [m]")
#     plot!(tref, vals_act, label=["x_act" "y_act" "z_act"])
    
    
end

