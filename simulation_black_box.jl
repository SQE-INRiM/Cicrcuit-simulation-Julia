#-------------------------------------SIMULATION AT LOW PUMP POWER-------------------------------------------

function simulate_low_pump_power(params, sim_vars, circuit, circuitdefs)

    
    # Initialize arrays for simulation output
    outvalsS21phidcSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:phidcSweep]))
    outvalsS12phidcSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:phidcSweep]))
    outvalsS11phidcSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:phidcSweep]))
    outvalsS22phidcSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:phidcSweep]))
    outvalsS21PhasephidcSweep = zeros(Float64, length(sim_vars[:ws]), length(sim_vars[:phidcSweep]))
    
    println("First simulation")
    
    # Run the simulation at low pump power
    @time for (k, phidcLocal) in enumerate(sim_vars[:phidcSweep])
        sources = [
            (mode=(0,), port=3, current=phidcLocal * (2 * 280 * 1e-6)),
            (mode=(1,), port=1, current=sim_vars[:Ip]),
        ]
        sol = hbsolve(
            sim_vars[:ws], sim_vars[:wp], sources, (1,), (1,),
            circuit, circuitdefs;
            dc=true, threewavemixing=true, fourwavemixing=true, iterations=sim_vars[:Niterations]
        )
        outvalsS21phidcSweep[:, k] = sol.linearized.S((0,), 2, (0,), 1, :)
        outvalsS12phidcSweep[:, k] = sol.linearized.S((0,), 1, (0,), 2, :)
        outvalsS11phidcSweep[:, k] = sol.linearized.S((0,), 1, (0,), 1, :)
        outvalsS22phidcSweep[:, k] = sol.linearized.S((0,), 2, (0,), 2, :)
        outvalsS21PhasephidcSweep[:, k] = unwrap(angle.(sol.linearized.S((0,), 2, (0,), 1, :)))
    end

    return outvalsS21phidcSweep, outvalsS12phidcSweep, outvalsS11phidcSweep, outvalsS22phidcSweep, outvalsS21PhasephidcSweep

end 


function plot_low_pump_power(params, sim_vars, circuit, circuitdefs)

    outvalsS21phidcSweep, outvalsS12phidcSweep, outvalsS11phidcSweep, outvalsS22phidcSweep, outvalsS21PhasephidcSweep = simulate_low_pump_power(params, sim_vars, circuit, circuitdefs)

    # Generate plots-----------------------------------------------------------------------

    p1 = plot(
        sim_vars[:phidcSweep],
        sim_vars[:ws] / (2 * pi * 1e9),
        10 * log10.(abs2.(outvalsS11phidcSweep)),
        seriestype=:heatmap,
        c=cgrad(:viridis, rev=true),
        clim=(-30, 0),
        xlabel=L"\phi / \phi_{0}",
        ylabel=L"f / GHz",
        title=L"S_{11} / dB",
        legend=false,
        colorbar=true
    )
    

    p2 = plot(
        sim_vars[:phidcSweep],
        sim_vars[:ws] / (2 * pi * 1e9),
        10 * log10.(abs2.(outvalsS21phidcSweep)),
        seriestype=:heatmap,
        c=:viridis,
        clim=(-5, 0),
        xlabel=L"\phi / \phi_{0}",
        ylabel=L"f / GHz",
        title=L"S_{21} / dB",
        legend=false,
        colorbar=true
    )

    p3 = plot(
        sim_vars[:phidcSweep],
        sim_vars[:ws] / (2 * pi * 1e9),
        20 * log10.(1 .- abs2.(outvalsS11phidcSweep) .- abs2.(outvalsS21phidcSweep)),
        seriestype=:heatmap,
        c=:viridis,
        clim=(-40, 0),
        xlabel=L"\phi / \phi_{0}",
        ylabel=L"f / GHz",
        title="Internal Losses / dB",
        legend=false,
        colorbar=true
    )


    phidcIndex = findall(x -> x == sim_vars[:phidc], sim_vars[:phidcSweep])

    p4 = plot(
        sim_vars[:ws] / (2 * pi * 1e9),
        -outvalsS21PhasephidcSweep[:, phidcIndex] / params[:N],
        xlabel=L"f / GHz",
        ylabel=L"k / rad \cdot cells^{-1}",
        title="Dispersion relation",
        #ylim=(0.0, 1.5),
        legend=true,
        colorbar=true,
        label="",
        framestyle=:box
    )

    vline!(p4, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
    vline!(p4, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
    vline!(p4, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
    vline!(p4, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")

    wp = 2*pi* round(sim_vars[:fp], digits=-8)                                    #Put the nearest wp value included inside ws. The reason of this line is because for computational reason the wp cannot be a value of ws.

    wpIndex = findall(x -> x == wp, sim_vars[:ws])
    wphalfIndex = findall(x -> x == wp/2, sim_vars[:ws])
    
    wp=sim_vars[:ws][wpIndex]
    wphalf=sim_vars[:ws][wphalfIndex]

    #linear relation
    y1=-outvalsS21PhasephidcSweep[10, phidcIndex] / params[:N]
    y2=-outvalsS21PhasephidcSweep[1, phidcIndex] / params[:N]
    x1=sim_vars[:ws][10]
    x2=sim_vars[:ws][1]

    m = (y1-y2)/(x1-x2)
    q = y2-m*x2
    
    plot!(p4, sim_vars[:ws] / (2 * pi * 1e9), m .* sim_vars[:ws] .+ q, label="linear relation", color=:orange)
    
    #line passing throug wp
    y1=-outvalsS21PhasephidcSweep[wpIndex[1], phidcIndex] / params[:N]
    y2=-outvalsS21PhasephidcSweep[1, phidcIndex] / params[:N]
    x1=sim_vars[:ws][wpIndex[1]]
    x2=sim_vars[:ws][1]

    m_p = (y1-y2)/(x1-x2)
    q_p = y2-m_p*x2
    
    plot!(p4, sim_vars[:ws] / (2 * pi * 1e9), m_p .* sim_vars[:ws] .+ q_p, label="wp line", color=:darkblue)
    

    #line passing throug wp/2 
    y1=-outvalsS21PhasephidcSweep[wphalfIndex[1], phidcIndex] / params[:N]
    y2=-outvalsS21PhasephidcSweep[1, phidcIndex] / params[:N]
    x1=sim_vars[:ws][wphalfIndex[1]]
    x2=sim_vars[:ws][1]

    m_phalf=(y1-y2)/(x1-x2)
    q_phalf = y2-m_phalf*x2
    
    plot!(p4, sim_vars[:ws] / (2 * pi * 1e9), m_phalf .* sim_vars[:ws] .+ q_phalf, label="wp/2 line", color=:darkred)

    alpha_wphalf=atan(m_phalf[1])
    alpha_wp=atan(m_p[1])
    alpha_lin=atan(m[1])


    return  p1, p2, p3, p4          

end


function calculation_low_pump_power(params, sim_vars, circuit, circuitdefs)

    outvalsS21phidcSweep, outvalsS12phidcSweep, outvalsS11phidcSweep, outvalsS22phidcSweep, outvalsS21PhasephidcSweep = simulate_low_pump_power(params, sim_vars, circuit, circuitdefs)

    wp = 2*pi* round(sim_vars[:fp], digits=-8)                                    #Put the nearest wp value included inside ws. The reason of this line is because for computational reason the wp cannot be a value of ws.

    phidcIndex = findall(x -> x == sim_vars[:phidc], sim_vars[:phidcSweep])
    wpIndex = findall(x -> x == wp, sim_vars[:ws])
    wphalfIndex = findall(x -> x == wp/2, sim_vars[:ws])
    
    wp=sim_vars[:ws][wpIndex]
    wphalf=sim_vars[:ws][wphalfIndex]

    #linear relation
    y1=-outvalsS21PhasephidcSweep[10, phidcIndex] / params[:N]
    y2=-outvalsS21PhasephidcSweep[1, phidcIndex] / params[:N]
    x1=sim_vars[:ws][10]
    x2=sim_vars[:ws][1]

    m = (y1-y2)/(x1-x2)
    q = y2-m*x2
    
    
    #line passing throug wp
    y1=-outvalsS21PhasephidcSweep[wpIndex[1], phidcIndex] / params[:N]
    y2=-outvalsS21PhasephidcSweep[1, phidcIndex] / params[:N]
    x1=sim_vars[:ws][wpIndex[1]]
    x2=sim_vars[:ws][1]

    m_p = (y1-y2)/(x1-x2)
    q_p = y2-m_p*x2
    

    #line passing throug wp/2 
    y1=-outvalsS21PhasephidcSweep[wphalfIndex[1], phidcIndex] / params[:N]
    y2=-outvalsS21PhasephidcSweep[1, phidcIndex] / params[:N]
    x1=sim_vars[:ws][wphalfIndex[1]]
    x2=sim_vars[:ws][1]

    m_phalf=(y1-y2)/(x1-x2)
    q_phalf = y2-m_phalf*x2

    alpha_wphalf=atan(m_phalf[1])
    alpha_wp=atan(m_p[1])
    alpha_lin=atan(m[1])

    return alpha_wphalf, alpha_wp, alpha_lin

end

#------------------------------------------SIMULATION AT FIXED FLUX---------------------------------------------------


function simulate_at_fixed_flux(sim_vars, circuit, circuitdefs)
    
    # Accessing ws inside the sim_vars dictionary
    outvalsS21IpSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:IpSweep]))
    outvalsS12IpSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:IpSweep]))
    outvalsS11IpSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:IpSweep]))
    outvalsS22IpSweep = zeros(Complex{Float64}, length(sim_vars[:ws]), length(sim_vars[:IpSweep]))


    println("Second simulation")

    @time for (k, IpSweepLocal) in enumerate(sim_vars[:IpSweep])
        sources = [
            (mode=(0,), port=3, current=sim_vars[:phidc] * (2 * 280 * 1e-6)),
            (mode=(1,), port=1, current=IpSweepLocal),
        ]
        
        sol = hbsolve(sim_vars[:ws], sim_vars[:wp], sources, sim_vars[:Nmodulationharmonics], sim_vars[:Npumpharmonics],
                    circuit, circuitdefs; dc=true, threewavemixing=true, fourwavemixing=true, iterations=sim_vars[:Niterations])
        outvalsS21IpSweep[:, k] = sol.linearized.S((0,), 2, (0,), 1, :)
        outvalsS12IpSweep[:, k] = sol.linearized.S((0,), 1, (0,), 2, :)
        outvalsS11IpSweep[:, k] = sol.linearized.S((0,), 1, (0,), 1, :)
        outvalsS22IpSweep[:, k] = sol.linearized.S((0,), 2, (0,), 2, :)
    end



    # Generate plots-----------------------------------------------------------------------------------------------
    
    cmap = cgrad([:darkblue, :white, :darkred])

    p1p = plot(
        sim_vars[:IpSweep] / 1e-6,
        sim_vars[:ws] / (2 * pi * 1e9),
        10 * log10.(abs2.(outvalsS11IpSweep)),
        seriestype=:heatmap,
        c = cmap,
        clim=(-20, 20),
        xlabel=L"I_{p} / \mu A",
        ylabel=L"f / GHz",
        legend=false,
        colorbar=true
    )

    p2p = plot(
        sim_vars[:IpSweep] / 1e-6,
        sim_vars[:ws] / (2 * pi * 1e9),
        10 * log10.(abs2.(outvalsS21IpSweep)),
        seriestype=:heatmap,
        c = cmap,
        clim=(-20, 20),
        xlabel=L"I_{p} / \mu A",
        ylabel=L"f / GHz",
        legend=false,
        colorbar=true
    )

    IpIndex = findall(x -> x == sim_vars[:IpGain], sim_vars[:IpSweep])

    p5 = plot(
        sim_vars[:ws] / (2 * pi * 1e9),
        10 * log10.(abs2.(outvalsS21IpSweep)[:, IpIndex]),
        xlabel=L"f / GHz",
        ylabel=L"S_{21} / dB",
        title="Gain Profile",
        legend=false,
        colorbar=true,
        framestyle=:box
    )

    return p1p, p2p, p5

end



#----------------------------------------------COSMESI FINALE-----------------------------------------------------------

function final_report(params, sim_vars, fixed_params, p1, p2, p3, p4, p1p, p2p, p5)

    
    # Plot customization
    vline!(p1, [sim_vars[:phidc]], width=2, style=:dash, color=:black)
    vline!(p2, [sim_vars[:phidc]], width=2, style=:dash, color=:black)
    vline!(p3, [sim_vars[:phidc]], width=2, style=:dash, color=:black)

    hline!(p1, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black)
    hline!(p2, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black)
    hline!(p3, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black)
    hline!(p1p, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black)
    hline!(p2p, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black)
    

    vline!(p4, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
    vline!(p4, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
    vline!(p4, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
    vline!(p4, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")


    vline!(p1p, [sim_vars[:IpGain] / 1e-6], width=2, color=:black)
    vline!(p2p, [sim_vars[:IpGain] / 1e-6], width=2, color=:black)

    hline!(p1p, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray)
    hline!(p2p, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray)
    hline!(p1p, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray)
    hline!(p2p, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray)
    hline!(p1p, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray)
    hline!(p2p, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray)

    vline!(p5, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black)

    value_fp = @sprintf("%.2f", sim_vars[:wp][1] / (2 * pi * 1e9))
    value_IpGain = @sprintf("%.4f", sim_vars[:IpGain] / (1e-6))
    value_Supercells = @sprintf("%.0f", params[:N] / params[:loadingpitch])
    value_criticalCurrentDensity = @sprintf("%.3f", params[:criticalCurrentDensity])
    value_CgDensity = @sprintf("%.1f", params[:CgDensity] / 1e-15)
    value_Cg = @sprintf("%.2f", params[:CgAreaUNLoaded] * params[:CgDensity] / 1e-15)


    empty_plot = plot([], legend=false, grid=false, framestyle=:none)  # Create an empty plot for the text

    
    annotate!(empty_plot, 0.0, 0.5, text("""
    Design parameters:
       loadingpitch = $(params[:loadingpitch]),
       A_small = $(params[:smallJunctionArea]) μm²,
       alphaSNAIL = $(params[:alphaSNAIL]),
       LloadingCell = $(params[:LloadingCell]),
       CgloadingCell = $(params[:CgloadingCell]),
       CgAreaUNLOADED = $(params[:CgAreaUNLoaded])
    
    Fabrication parameters:
       Jc = $(value_criticalCurrentDensity) μA/μm²,
       Cj density = $(fixed_params[:JosephsonCapacitanceDensity]) fF/μm²,
       CgDielectricThickness = $(params[:CgDielectricThichness]),
       Cg density = $(value_CgDensity) fF/μm²,
       Cg = $(value_Cg) fF
    
    Simulation parameters:
       Number of Cells = $(params[:N]),
       Number of Supercells = $(value_Supercells),
       fp = $(value_fp) GHz,
       Ip = $(value_IpGain) μA,
       Phi_dc = $(sim_vars[:phidc]),
    
       NPumpHarm = $(sim_vars[:Npumpharmonics][1]),
       NModHarm = $(sim_vars[:Nmodulationharmonics][1]),
       NIter = $(sim_vars[:Niterations])
    """, 12, :black, halign=:left))
    


    p = plot(p1, p2, p1p, p2p, p3, p4, p5, empty_plot, layout=(4, 2), margin=10Plots.mm, size=(1500, 2000))
    
    return p


end


#-----------------------------GENERAL FUNCTION FOR SIMULATION AND PLOT--------------------------------------

function simulate_and_plot(params_temp, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)

    p1,p2,p3,p4 = plot_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp)

    p1p,p2p,p5 = simulate_at_fixed_flux(sim_vars, circuit_temp, circuitdefs_temp)

    p_temp = final_report(params_temp, sim_vars, fixed_params, p1, p2, p3, p4, p1p, p2p, p5)

    return p_temp

end