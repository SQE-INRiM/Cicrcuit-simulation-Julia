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

using Surrogates
using LinearAlgebra
using QuasiMonteCarlo

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
    :JosephsonCapacitanceDensity => 10, # 45 Unit value
    :CgDielectricK => 9.6,             # Dielectric constant
    :tandeltaCgDielectric => 2.1e-3,   # Loss tangent of dielectric
)

#SIMULATION PARAMETERS

sim_vars = Dict(
    
    :fs => (3:0.1:30.0) * 1e9,
    :fp =>  13.301 * 1e9,
    :Ip => 0.00001e-6,
    :IpGain => 0.25e-6,
    :IpSweep => (0.1:0.05:0.3)*1e-6, #0.8*1e-6
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


"""
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


"""


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
                    
                    
optimal_parmas = add_parameters(optimal_params)

#possible benchmark
optimal_params = Dict(:CgDensity => 6.488427480916029e-16, :smallJunctionArea => 0.5, :CgDielectricThichness => 131.0, :alphaSNAIL => 0.175, :CgloadingCell => 1.5, :criticalCurrentDensity => 0.4, :phidc => 0.37, :nMacrocells => 150.0, :LloadingCell => 1.25, :N => 300.0, :CgAreaUNLoaded => 200.0, :loadingpitch => 3.0)

params=optimal_params
println("Optimal Parameters: $optimal_params")
#println("Optimal Metric: $optimal_metric")


circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, optimal_params, fixed_params)

S21, _, S11, _, S21phase = simulate_low_pump_power(sim_vars, circuit_temp, circuitdefs_temp)

p1, _, _, _ = plot_low_pump_power(S21, S11, S21phase, optimal_params, sim_vars)



#-----------------------------------------------------------------------------------


phidcIndex = findall(x -> x == params[:phidc], sim_vars[:phidcSweep])

y=-S21phase[:, phidcIndex] / params[:N]
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

using SavitzkyGolay

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


y_der_sg = savitzky_golay(y[:,1], window_size, poly_order, deriv=1)
#y_der_sg_rate = savitzky_golay(y[:,1], window_size, poly_order, deriv=1, rate=200/(15-(-5)))

p_sg_der = plot(
    x,
    [y_der_sg.y],
    xlabel=L"f / GHz",
    ylabel=L"k / rad \cdot cells^{-1}",
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

y_der2_sg = savitzky_golay(y[:,1], window_size, poly_order, deriv=2)
#y_der_sg_rate = savitzky_golay(y[:,1], window_size, poly_order, deriv=1, rate=200/(15-(-5)))

p_sg_der2 = plot(
    x,
    [y_der2_sg.y],
    xlabel=L"f / GHz",
    ylabel=L"k / rad \cdot cells^{-1}",
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



#---------------------------------------------------------------------------

p=plot(p1,p4, p_sm, p_sg, p_sg_der, p_sg_der2,  layout=(6,1), size=(1200, 1400))
display(p)

#maxS11 = maxS11val_BandFreq_FixFlux(S11, optimal_params, sim_vars)
#println("Maximum S11: ", maxS11)

p_temp = simulate_and_plot(optimal_params, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)
display(p_temp)

