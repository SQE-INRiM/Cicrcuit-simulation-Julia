#------------------------------------PACKAGES---------------------------------------------------------
using JosephsonCircuits
using Plots
using Random
using DSP
using Printf
using Dates
using LaTeXStrings
using ColorTypes
using Symbolics
using Dates
using Distributions
using NaNStatistics
using FindPeaks1D
using Surrogates
using LinearAlgebra
using QuasiMonteCarlo
using JLD2
using SavitzkyGolay

include("snail_circuit.jl")
include("simulation_black_box.jl")
include("utils.jl")



"""
#-----------------------------------LIST OF PARAMETERS--------------------------------------------------
#-------------------------------------------------------------------------------------------------------

#DESIGN PARAMETERS
loadingpitch = 2                #Tra 2 e 4 step 1!!                 #Indica pacchetto intero, c'è sempre 1 unloaded e n-1 loaded. Es: loadingpitch=3 --> 1 unloaded e 2 loaded
smallJunctionArea = 1           #Tra 2 e 5 (meglio basso)   
alphaSNAIL = 0.2                #Tra 0.1 e 0.25            
LloadingCell = 2                #Tra 1 e 6
CgloadingCell = 1               #Tra 1 e 6

#DESIGN PARAMETERS FIXED 
CgAreaUNLoaded = 150+20(smallJunctionArea/alphaSNAIL)                         #um^2 - Valore unitario,
nodePerCell = 4                            

#------------------------------------------------------------------------------------------------------

#FABRICATION PARAMETERS
JJSmallStd = 0.05               #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread
JJBigStd = 0.05                 #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread

criticalCurrentDensity = 1      #Tra 0.2 e 1         
CgDielectricThichness = 80      #Tra 10 e 70       


#FABRICATION PARAMETERS FIXED
JosephsonCapacitanceDensity = 45            

CgDielectricK = 9.6                         
tandeltaCgDielectric = 2.1e-3               

CgDensity = (CgDielectricK*8.854e-12)/((1e12)*CgDielectricThichness*(1e-9)) #F/um^2) 


#-------------------------------------------------------------------------------------------------------


#SIMULATION PARAMETERS

N=3*loadingpitch                           #number of simulated cells - 440 in https://doi.org/10.1063/5.0127690 (Roudsari - Chalmers paper)
ws=2*pi*(3.0:0.1:30.0)*1e9                  #Hz - signal range (start:step:stop) 
wp=(2*pi*8.001*1e9,)                         #Hz - pump frequency
Ip=0.00001e-6                                #A  - pump current (used for phidc Sweep simulation)
IpGain=0.8e-6
IpSweep=(0.8:0.1:0.8)*1e-6                  #A  - pump current (start:step:stop) NO SWEEP per ora

phidc=0.38                                 #flux induced by the flux line (used for Ip Sweep simulation)
phidcSweep=(0:0.01:0.5)                      #flux sweep induced by the flux line (start:step:stop)

Npumpharmonics = (8,)                        #used in IpSweep (1 in phidcSweep) - 8 in arXiv:2408.17293v1
Nmodulationharmonics = (4,)                  #used in IpSweep (1 in phidcSweep) - 4 in arXiv:2408.17293v1

Niterations=500                             #simulation max iterations (500 fine, 1000 ultrafine)


#----------------------------------------------------------------------------------


#SIMULATION CIRCUIT PARAMETERS FIXED

Rleft = 50.0,                                                          #Ohm
Rright = 50.0,                                                         #Ohm
Ladd = 70.0e-15,                                                       #Henry - loop inductance - Valore preso da TableI di arXiv:2408.17293v1
Lg = 20.0e-9,                                                          #Henry - geometrical inductance - Valore preso da TableI di arXiv:2408.17293v1
Lf = 190.0e-12,                                                        #Henry - Flux line inductors - Valore preso da TableI di arXiv:2408.17293v1
Cf = 0.076e-12,                                                        #Farad - Flux line capacitors - Valore preso da TableI di arXiv:2408.17293v1
M = 0.999,                                                             # the inverse inductance matrix for K=1.0 diverges, so set K<1.0

"""


#------------------------------------SIMULATION---------------------------------------------------------
#-------------------------------------------------------------------------------------------------------


#FIXED PARAMS

fixed_params = Dict(

    # Design parameters (fixed)
    :nodePerCell => 4,                 # Nodes per cell

    # Fabrication parameters (fixed)
    :JosephsonCapacitanceDensity => 10, #Unit value
    :CgDielectricK => 9.6,             # Dielectric constant
    :tandeltaCgDielectric => 2.1e-3,   # Loss tangent of dielectric
)

#SIMULATION PARAMETERS

sim_vars = Dict(
    
    :fs => (3:0.1:30.0) * 1e9,
    :fp =>  13.301 * 1e9,
    :Ip => 0.00001e-6,
    :IpGain => 0.25e-6,
    :IpSweep => (0:0.2:1)*1e-6, #(0.1:0.1:0.3)*1e-6, # #0.8*1e-6
    :phidcSweep => (0:0.01:0.5),
    :Npumpharmonics => (8,),
    :Nmodulationharmonics => (4,),
    :Niterations => 500
)

sim_vars[:ws] = 2 * pi * sim_vars[:fs]
sim_vars[:wp] = (2 * pi * sim_vars[:fp],)


#--------------------------------------------------------------------------------


#PARAMETER SCAN

JJSmallStd = 0.0               #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread
JJBigStd = 0.0                 #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread



# Define the parameter arrays directly in the dictionary
sim_params_space = Dict(
    :loadingpitch => [2],
    :nMacrocells => [20],                        
    :smallJunctionArea => collect(2:0.5:5),                        
    :alphaSNAIL => collect(0.1:0.05:0.25),                         
    :LloadingCell => [1, 2, 3],                           
    :CgloadingCell => [1, 2],                             
    :criticalCurrentDensity => [0.2, 0.5],                
    :CgDielectricThichness => [10, 20, 40, 70]                    
)





#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------


result=((1.0003413536085437, 2.0, 2.0006331298076976, 0.4646196724801679, 20.0, 3.8990685664598863, 69.60013491788783, 0.17201051536816814), 2.873839683907934)
keys_list = collect(keys(sim_params_space))  # Extract parameter names
optimal_vec = result[1]                # Optimized vector
optimal_metric = result[2]             # Optimal metric value

optimal_params = vector_to_param(optimal_vec, keys_list)





println("-----------------------------------------------------")
println("-----------------------------------------------------")
println("------------------SIMULATION RESULTS-----------------")
println("-----------------------------------------------------")
println("-----------------------------------------------------")


optimal_params = Dict(:smallJunctionArea => 0.7, 
                    :CgDielectricThichness => 90, 
                    :alphaSNAIL => 0.16, 
                    :CgloadingCell => 1, 
                    :criticalCurrentDensity => 0.4, 
                    :phidc => 0.38, 
                    :nMacrocells => 50, 
                    :LloadingCell => 1.5,
                    :loadingpitch => 3.0)
                    
                    
optimal_params = add_parameters(optimal_params)


#optimal_params=Dict(:CgDensity => 1.2635933809032367e-15, :smallJunctionArea => 2.1394323878236645, :CgDielectricThichness => 67.26720896499296, :alphaSNAIL => 0.15, :CgloadingCell => 1.5456783925252877, :criticalCurrentDensity => 0.4, :phidc => 0.36, :nMacrocells => 50.0, :LloadingCell => 1.5250587114918863, :N => 150.0, :CgAreaUNLoaded => 200.0, :loadingpitch => 3.0)
#possible benchmark
optimal_params = Dict(:CgDensity => 1.216548036096446e-15, :smallJunctionArea => 2.039918401851307, :CgDielectricThichness => 69.86851112984859, :alphaSNAIL => 0.15241400139422023, :CgloadingCell => 1.503267171440983, :criticalCurrentDensity => 0.4, :phidc => 0.36, :nMacrocells => 50.0, :LloadingCell => 1.5032002776845013, :N => 150.0, :CgAreaUNLoaded => 200.0, :loadingpitch => 3.0)
optimal_params=optimal_params
println("Optimal Parameters: $optimal_params")
#println("Optimal Metric: $optimal_metric")




circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, optimal_params, fixed_params)

#S21, _, S11, _, S21phase = simulate_low_pump_power(sim_vars, circuit_temp, circuitdefs_temp)

#p1, _, _, _ = plot_low_pump_power(S21, S11, S21phase, optimal_params, sim_vars)



#-----------------------------------------------------------------------------------


phidcIndex = findall(x -> x == optimal_params[:phidc], sim_vars[:phidcSweep])

y=-S21phase[:, phidcIndex] / optimal_params[:N]
x= sim_vars[:ws] / (2 * pi * 1e9)

p4 = plot(
    x,
    y,
    xlabel=L"f / GHz",
    ylabel=L"k / rad \cdot cells^{-1}",
    title="Dispersion relation",
    ylim=(0.0, 1.5),
    legend=true,
    colorbar=true,
    label="",
    framestyle=:box
)

vline!(p4, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
vline!(p4, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p4, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p4, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")


wp = 2*pi* round(sim_vars[:fp], digits=-8)                              #Put the nearest wp value included inside ws. The reason of this line is because for computational reason the wp cannot be a value of ws.
wphalf = 2*pi* round(sim_vars[:fp]/2, digits=-8)  


phidcIndex = findall(x -> x == optimal_params[:phidc], sim_vars[:phidcSweep])
wpIndex = findall(x -> x == wp, sim_vars[:ws])
wphalfIndex = findall(x -> x == wphalf, sim_vars[:ws])


#----------------------------------------------------------------------

#ks + ki = kp

f_p=x[wpIndex[1]]
k_p=y[wpIndex[1]]
println("f_p ", f_p)
println("k_p ", k_p)

f_phalf=x[wphalfIndex[1]]
k_phalf=y[wphalfIndex[1]]
println("f_phalf ", f_phalf)
println("k_phalf ", k_phalf)


fp_band_dx = round((sim_vars[:fp]/2)+1e9, digits=-8)  
fp_band_sx = round((sim_vars[:fp]/2)-1e9, digits=-8)
k_band_dx_index = findall(x -> x == fp_band_dx, sim_vars[:fs])
k_band_sx_index = findall(x -> x == fp_band_sx, sim_vars[:fs])
fp_band = x[k_band_sx_index[1]:k_band_dx_index[1]]  
kp_band_dx = y[k_band_dx_index[1]]
kp_band_sx = y[k_band_sx_index[1]]
kp_band = y[k_band_sx_index[1]:k_band_dx_index[1]]
println("f of f_band ", fp_band_sx, " ", fp_band_dx)
println("f band ", fp_band)
println("k of f_band ", kp_band_sx, " ", kp_band_dx)
println("k band ", kp_band)

function compute_deltaK(kp_band, k_p) 
    sums = 0
    n = length(kp_band)

    for idx in 1:div(n, 2)
        # Add symmetric pairs
        sums += kp_band[idx] + kp_band[end - idx + 1]
    end

    # Handle the middle element if the length of kp_band is odd
    if isodd(n)
        sums += kp_band[div(n, 2) + 1]
    end

    println(sums/length(kp_band))
    deltaK = sums/length(kp_band) - k_p
    println("deltaK ", deltaK)
    
    return abs(deltaK) 

end    

deltaKK = compute_deltaK(kp_band, k_p)
println("deltaKK ", deltaKK)

#deltaK = k_s + k_i - k_p





#------------------------------------------------------------------------



y_stopband = y[wpIndex[1]-10:wpIndex[1]-2]
x_stopband = sim_vars[:ws][wpIndex[1]-10:wpIndex[1]-2]

n = length(x_stopband)
X = [ones(n) x_stopband] # Add a column of ones for the intercept

beta = X \ y_stopband # Solves for [intercept, slope]

q_stopband = beta[1]
m_stopband = beta[2]

fitted_y = X * beta

scatter!(p4, x_stopband/(2 * pi * 1e9), y_stopband, label="Data Points")
plot!(p4, x_stopband/(2 * pi * 1e9), fitted_y, label="Linear Fit", lw=2)
plot!(p4, sim_vars[:ws] / (2 * pi * 1e9), m_stopband .* sim_vars[:ws] .+ q_stopband, label="stopband line", color=:darkred)




#-------------------------------------------------------------



window_size = 7 # Number of points for smoothing
smoothed_y = movmean(y, window_size)


p_sm = plot(
    x,
    smoothed_y,
    xlabel=L"f / GHz",
    ylabel=L"k / rad \cdot cells^{-1}",
    title="DR with move mean",
    ylim=(0.0, 1.5),
    legend=true,
    colorbar=true,
    label="",
    framestyle=:box
)

vline!(p_sm, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
vline!(p_sm, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sm, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sm, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")

y_stopband = smoothed_y[wpIndex[1]-10:wpIndex[1]-2]
x_stopband = sim_vars[:ws][wpIndex[1]-10:wpIndex[1]-2]

n = length(x_stopband)
X = [ones(n) x_stopband] # Add a column of ones for the intercept

beta = X \ y_stopband # Solves for [intercept, slope]

q_stopband = beta[1]
m_stopband = beta[2]

fitted_y = X * beta

scatter!(p_sm, x_stopband/(2 * pi * 1e9), y_stopband, label="Data Points")
plot!(p_sm, x_stopband/(2 * pi * 1e9), fitted_y, label="Linear Fit", lw=2)
plot!(p_sm, sim_vars[:ws] / (2 * pi * 1e9), m_stopband .* sim_vars[:ws] .+ q_stopband, label="stopband line", color=:darkred)

#---------------------------------------------------------------------------


window_size = 11 # Must be odd
poly_order = 3
y_sg = savitzky_golay(y[:,1], window_size, poly_order)    

p_sg = plot(
    x,
    [y_sg.y],
    xlabel=L"f / GHz",
    ylabel=L"k / rad \cdot cells^{-1}",
    title="DR with SavitzkyGolay",
    ylim=(0.0, 1.5),
    legend=true,
    colorbar=true,
    label="",
    framestyle=:box
)


vline!(p_sg, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
vline!(p_sg, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")

y_stopband = y_sg.y[wpIndex[1]-10:wpIndex[1]-2]
x_stopband = sim_vars[:ws][wpIndex[1]-10:wpIndex[1]-2]

n = length(x_stopband)
X = [ones(n) x_stopband] # Add a column of ones for the intercept

beta = X \ y_stopband # Solves for [intercept, slope]

q_stopband = beta[1]
m_stopband = beta[2]

fitted_y = X * beta

scatter!(p_sg, x_stopband/(2 * pi * 1e9), y_stopband, label="Data Points")
plot!(p_sg, x_stopband/(2 * pi * 1e9), fitted_y, label="Linear Fit", lw=2)
plot!(p_sg, sim_vars[:ws] / (2 * pi * 1e9), m_stopband .* sim_vars[:ws] .+ q_stopband, label="stopband line", color=:darkred)

#---------------------------------------------------------------------------

window_size = 21 # Must be odd
poly_order = 3

y_der_sg = savitzky_golay(y[:,1], window_size, poly_order, deriv=1, rate=100)
#y_der_sg_rate = savitzky_golay(y[:,1], window_size, poly_order, deriv=1, rate=200/(15-(-5)))

p_sg_der = plot(
    x,
    [y_der_sg.y],
    xlabel=L"f / GHz",
    ylabel=L"a. u.",
    title="First derivative DR with SavitzkyGolay",
    #ylim=(0.0, 1.5),
    legend=true,
    colorbar=true,
    label="",
    framestyle=:box
)


vline!(p_sg_der, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
vline!(p_sg_der, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg_der, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg_der, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")


#--------------------------------------------------------------------------

window_size = 11 # Must be odd
poly_order = 3

y_der2_sg = savitzky_golay(y[:,1], window_size, poly_order, deriv=2, rate=100)
#y_der2_sg_rate = savitzky_golay(y[:,1], window_size, poly_order, deriv=1, rate=200/(15-(-5)))

p_sg_der2 = plot(
    x,
    [y_der2_sg.y],
    xlabel=L"f / GHz",
    ylabel=L"a. u.",
    title="Second derivative DR with SavitzkyGolay",
    #ylim=(0.0, 1.5),
    legend=true,
    colorbar=true,
    label="",
    framestyle=:box
)


vline!(p_sg_der2, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
vline!(p_sg_der2, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg_der2, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg_der2, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")


y = y_der2_sg.y

"""
id_lb_range1 = findall(x -> x == 12, x)[1]
id_ub_range1 = findall(x -> x == 13, x)[1]

max_y, index = findmax(y[id_lb_range1:id_ub_range1])
max_x = x[id_lb_range1:id_ub_range1][index]
#scatter!(p_sg_der2, [max_x], [max_y], label="Max")


id_lb_range2 = findall(x -> x == 11, x)[1]
id_ub_range2 = findall(x -> x == 13, x)[1]

min_y, index = findmin(y[id_lb_range2:id_ub_range2])
min_x = x[id_lb_range2:id_ub_range2][index]
#scatter!(p_sg_der2, [min_x], [min_y], label="Min")
"""
#Finding min --> ci interessa il max ma il primo min è piu sensibile

pkindices, properties = findpeaks1d(-y; 
    height=5,             # Minimum peak height
    prominence=0.01,        # Minimum prominence of peaks
    width=1.0,             # Minimum width of peaks
    relheight=0.5           # Relative height to determine peak edges
)

x_mins = x[pkindices]
y_mins = y[pkindices]
scatter!(p_sg_der2, x_mins, y_mins, color="blue", markersize=2, label="Min")



#println("xmins: ", x_mins)
#first_x_min = isempty(x_mins) ? 0 : x_mins[1]
#println("first_x_min ", first_x_min)




# Finding max

pkindices, properties = findpeaks1d(y; 
    height=1,             # Minimum peak height
    prominence=0.01,        # Minimum prominence of peaks
    width=1.0,             # Minimum width of peaks
    relheight=0.5           # Relative height to determine peak edges
)

x_maxs = x[pkindices]
y_maxs = y[pkindices]

scatter!(p_sg_der2, x_maxs, y_maxs, color="red", markersize=2, label="Max")



function find_first_peak(x_maxs, y_maxs, first_x_min)
    if first_x_min == 0
        return 0, 0
    end
    for (x_max, y_max) in zip(x_maxs, y_maxs)
        if x_max > first_x_min 
            return x_max, y_max
        end
    end
    return 0,0  # Return nothing if no peak is found
end

x_peak_sb, y_peak_sb = find_first_peak(x_maxs, y_maxs, first_x_min)
#println(x_peak_sb)
#println( y_peak_sb)

scatter!(p_sg_der2, [x_peak_sb], [y_peak_sb], color="green", markersize=4, label="Stopband peak")





#---------------------------------------------------------------------------

























"""
window_size = 11 # Must be odd
poly_order = 5

y_der3_sg = savitzky_golay(y[:,1], window_size, poly_order, deriv=3)

p_sg_der3 = plot(
    x,
    [y_der3_sg.y],
    xlabel=L"f / GHz",
    ylabel=L"k''' / a. u.",
    title="Third derivative DR with SavitzkyGolay",
    #ylim=(0.0, 1.5),
    legend=true,
    colorbar=true,
    label="",
    framestyle=:box
)


vline!(p_sg_der3, [sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:black, label="")
vline!(p_sg_der3, [(1 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg_der3, [(3 / 2) * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, style=:dash, color=:gray, label="")
vline!(p_sg_der3, [2 * sim_vars[:wp][1] / (2 * pi * 1e9)], width=2, color=:gray,label="")

dy=y_der3_sg.y

# Find zeros of the derivative (smoothed version)
function find_der_max(x, y, dy; neighborhood=1, prominence_threshold=1e-5)
    extrema = []
    extrema_dx = []
    n = length(dy)
    for i in (neighborhood + 1):(n - neighborhood)
        
        if dy[i-1] > 0 && dy[i] < 0
            # Check if it's a significant maximum
            is_prominent = all(y[i] >= y[j] + prominence_threshold for j in (i-neighborhood):(i+neighborhood) if j != i)
            
            if is_prominent
                push!(extrema, (x[i], y[i]))  # Local maximum
                push!(extrema_dx, (x[i], dy[i])) 
                
            end
        end
    end
    return extrema, extrema_dx
end

#extrema, extrema_dx = find_der_max(x, y, dy)
maxima = findlocalmaxima(imfilter(y, KernelFactors.gaussian((3))))
maxima = [getindex(max, 1) for max in maxima]

# Output extrema
println(extrema)

# Optional: Plot the results
scatter!(p_sg_der2, [e[1] for e in extrema], [e[2] for e in extrema], label="Extrema", color=:red)
#hline!(p_sg_der3, [0], width=1, color=:black, label="")
#scatter!(p_sg_der3, [e[1] for e in extrema_dx], [e[2] for e in extrema_dx], label="Extrema", color=:blue)
"""

#----------------------------------------------------------------------------

p=plot(p4, p_sg, p_sg_der, p_sg_der2, layout=(4,1), size=(1200, 1400))
display(p)

#maxS11 = maxS11val_BandFreq_FixFlux(S11, optimal_params, sim_vars)
#println("Maximum S11: ", maxS11)

#p_temp = simulate_and_plot(optimal_params, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)
#display(p_temp)
