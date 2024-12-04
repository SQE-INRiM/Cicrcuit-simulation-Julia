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
using Distributions
using Dates
using StatsBase
using SavitzkyGolay
using FindPeaks1D
using Surrogates
using LinearAlgebra
using QuasiMonteCarlo
using JLD2

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
    :JosephsonCapacitanceDensity => 45, # Unit value
    :CgDielectricK => 9.6,             # Dielectric constant
    :tandeltaCgDielectric => 2.1e-3,   # Loss tangent of dielectric
)

#SIMULATION PARAMETERS

sim_vars = Dict(
    
    :fs => (3:0.1:30.0) * 1e9,
    :fp =>  14.001*1e9, #13.301 * 1e9,
    :Ip =>  0.00001e-6,
    :IpSweep => (0:0.1:2.5)*1e-6, #0.8*1e-6 iniziale
    :phidcSweep => (0:0.01:0.5),
    :Npumpharmonics => (8,),
    :Nmodulationharmonics => (4,),
    :Niterations => 500
)

sim_vars[:ws] = 2 * pi * sim_vars[:fs]
sim_vars[:wp] = (2 * pi * sim_vars[:fp],)


#--------------------------------------------------------------------------------


#PARAMETER SCAN

JJSmallStd = 0.1        #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread
JJBigStd = 0.1               #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread

# Define the parameter arrays directly in the dictionary
"""
# Parameter space near the benchmark
sim_params_space = Dict(
    :loadingpitch => [3],
    :nMacrocells => [50],   #20                     
    :smallJunctionArea =>  collect(0.5:0.25:1.5), # [0.7]                        
    :alphaSNAIL => collect(0.1:0.025:0.2),  #0.16                  
    :LloadingCell => collect(1:0.25:2),   #[1.5], #                        
    :CgloadingCell => collect(0.5:0.25:1.5),  #[1], #                           
    :criticalCurrentDensity => [0.4], #collect(0.1:0.1:0.9),                
    :CgDielectricThichness => collect(50:1:140)     
 
)

"""

# Right parameter space

sim_params_space = Dict(
    :loadingpitch => [3],
    :nMacrocells => [50],               
    :smallJunctionArea =>  collect(0.5:0.5:5),                       
    :alphaSNAIL => collect(0.1:0.05:0.25),                   
    :LloadingCell => collect(1:1.5:4.5),                          
    :CgloadingCell => collect(1:1.5:4.5),                             
    :criticalCurrentDensity => [0.4],                 
    :CgDielectricThichness => collect(10:10:70)     
)

sim_params_space = Dict(
    :loadingpitch => [3],
    :nMacrocells => [50],               
    :smallJunctionArea =>  [2,3,4],                       
    :alphaSNAIL => collect(0.15:0.5:0.25),                   
    :LloadingCell => collect(1.5:1.5:4.5),                          
    :CgloadingCell => collect(1.5:1.5:4.5),                             
    :criticalCurrentDensity => [0.4],                 
    :CgDielectricThichness => collect(10:20:70)     
)

"""
# test par space
sim_params_space = Dict(
    :loadingpitch => [3],
    :nMacrocells => [50],               
    :smallJunctionArea =>  [2,3],                       
    :alphaSNAIL => [0.2],                   
    :LloadingCell => [3],                          
    :CgloadingCell => [3],                             
    :criticalCurrentDensity => [0.4],                 
    :CgDielectricThichness => [40]     
)

"""


#-----------------------------------------------------------------------------------
#-----------------------------STARTING SIMULATION-----------------------------------
#-----------------------------------------------------------------------------------

"""
First metric

#metric_angles = abs(delta_alpha_wp-delta_alpha_wphalf) * (1/((abs(delta_alpha_wp*delta_alpha_wphalf))^(1/2)))
metric_angles_stopband = (abs(alpha_stopband)*2.5)*1e12
println("   a. Angle contribute stopband: ", metric_angles_stopband)

metric_angles_nonlin = 1/(sqrt(abs(alpha_lin - alpha_wp))*1e6) # abs(alpha_nonlin - alpha_wp)*1e12
println("   b. Angles contribute non-linarity: ", metric_angles_nonlin)

metric_impedance = (20/abs(maxS11value))
println("   c. Impedance matching contibute: ", metric_impedance)

metric = metric_angles_stopband + metric_angles_nonlin + metric_impedance   
println("Final value: ", metric)

     #metric_angles = abs(delta_alpha_wp-delta_alpha_wphalf) * (1/((abs(delta_alpha_wp*delta_alpha_wphalf))^(1/2)))
  

    metric_stopband_position =2/(abs(min_y * max_y)) tra l'intervallo prima della fp

"""




# Define the cost function
function cost(vec) #vector passing
    
    println("-----------------------------------------------------")
    #println("Vector: ", vec)

    params_temp = vector_to_param(vec, keys(sim_params_space))
    #println("Dictionary: ", params_temp)

    params_temp = add_parameters(params_temp)
    println("Dictionary update: ", params_temp)
    #println("flux value: ", params_temp[:phidc], " for alpha: ", params_temp[:alphaSNAIL])

    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)
    println("Circuit created")

    _, _, S11, _, S21phase = simulate_low_pump_power(sim_vars, circuit_temp, circuitdefs_temp)
    println("S parameters calculated")

    maxS11value = maxS11val_BandFreq_FixFlux(S11, params_temp, sim_vars)
    println("Maximum value of S11: ", maxS11value)

    alpha_wphalf, alpha_wp, alpha_lin, alpha_stopband, alpha_nonlin  = calculation_metric_lines(S21phase, params_temp, sim_vars)
    println("Angles calculated")

    _, x_stopband_peak, x_pump = plot_derivative_low_pump(S21phase, sim_vars, params_temp)
    println("Peak calculated")

    #delta_alpha_wp=abs(alpha_wp - alpha_lin)
    #delta_alpha_wphalf=abs(alpha_wphalf - alpha_lin)

    metric_angles_stopband = (abs(alpha_stopband)*2.5)*1e12
    println("   a. Angle contribute stopband: ", metric_angles_stopband)

    metric_stopband_position = 3.5*abs(x_stopband_peak - x_pump)
    println("   b. Stopband position contribute: ", metric_stopband_position)

    metric_impedance = (20/abs(maxS11value))
    println("   c. Impedance matching contibute: ", metric_impedance)

    metric = metric_angles_stopband + metric_stopband_position + metric_impedance   
    println("Final value: ", metric)

    #p1,_,_,p4 = plot_low_pump_power(S21, S11, S21phase, params_temp, sim_vars)
    #p=plot(p1,p4,layout=(2,1), size=(600, 700))
    
    #display(p)

    println("-----------------------------------------------------")
  
    return metric
    #p_temp = simulate_and_plot(params_temp, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)
    #display(p_temp)

end



println("-----------------------------------------------------")

n_initial_points = find_n_initial_points(sim_params_space) 
n_maxiters = 10 #50
n_num_new_samples = 12 #80

time_estimated, finish_time = simulation_time_estimation(n_initial_points, n_maxiters, n_num_new_samples)

println("Total number of initial points: ", n_initial_points)
println("Number of max interactions: ", n_maxiters)
println("Number of samples for every interaction: ", n_num_new_samples)
println("Total number of point calculated in the simulation: ", n_maxiters*n_num_new_samples)

println("SIMULATION TIME ESTIMATED: ", time_estimated)

println("-----------------------------------------------------")
start_time = time()
dates_now = (Dates.now())
println("STARTING AT: ", dates_now)
println("ESTIMATED END AT: ", finish_time)
println("-----------------------------------------------------")

output = "SIMULATION IN PROGRESS... \n\nTIME ESTIMATED:  $(time_estimated)\n\n\n\nSTARTING AT: $(dates_now)\n\nESTIMATED END AT: $(finish_time)"
out_plot = plot([], legend=false, grid=false, framestyle=:none)
annotate!(out_plot, 0.5, 0.5, text(output, 12))  # Position the text at the center of the plot
display(out_plot)




#Return upper and lower boundaries
bounds = [extrema(v) for v in values(sim_params_space)]
println("Bounds: ", bounds)

lb = [b[1] for b in bounds]
ub = [b[2] for b in bounds]
println("Lower bounds:", lb)
println("Upper bounds:", ub)


# Create n random points

#n_initial_points = 2
#ini_points = generate_n_initial_random_points(n_initial_points, sim_params_space)
#println("Initial points: ", ini_points)


# Create all the points of the parameter space

initial_points = generate_all_initial_points(sim_params_space)
save_points_to_file(initial_points, "initial_points.jld2")
#println("Initial points: ", initial_points)

#for p in initial_points
#    println("Single point: ", p)
#end 

initial_values = [cost(p) for p in initial_points]
save_points_to_file(initial_values, "initial_values.jld2")
#println("Initial values: ", initial_values)


#initial_points = load_points_from_file("initial_points.jld2")
#println("Initial points: ", initial_points)

#initial_values = load_points_from_file("initial_values.jld2")
#println("Initial values: ", initial_values)


my_k_SRBFN = Kriging(initial_points, initial_values, lb, ub)
#println("my_k_SRBFN: ", my_k_SRBFN)





println("-----------------------------------------------------")
println("-----------------------------------------------------")
println("---STARTING GAUSSIAN PROCESS BAYESIAN OPTIMIZATION---")
println("-----------------------------------------------------")


# Perform Gaussian Process optimization

result = surrogate_optimize(
    cost,
    SRBF(),
    lb,
    ub,
    my_k_SRBFN,
    RandomSample(),
    maxiters = n_maxiters,                        # Number of interactions. Incresing maxiters: Leads to a longer optimization process with potentially better solutions but at the cost of more time.
    num_new_samples = n_num_new_samples     # Number of point generated for every single interaction. Incresing num_new_samples: Allows each iteration to consider a broader range of candidate points, 
                                            # improving the chance of finding a good solution early but also increasing the computational cost per iteration.
)


println("-----------------------------------------------------")
println("FINISHED AT: ", (Dates.now()))
println("Complete simulation time: ", round((time() - start_time) / 60), " minutes and ", Int(round((time()-start_time) % 60))," seconds.")

"""
Balancing these maxiters and num_new_samples is essential for efficient optimization. 
For example, a small maxiters and a large num_new_samples may work well for problems where sampling is expensive but model updates are less so.

If the problem has many variables or the surrogate model is very non-linear, you may want:
Larger num_new_samples: To explore the space better per iteration.
Smaller maxiters: To limit the total number of model updates.

If the problem is relatively simple or smooth:
Smaller num_new_samples: To reduce computational cost per iteration.
Larger maxiters: To allow for more refinement.
"""



keys_list = collect(keys(sim_params_space))  # Extract parameter names
optimal_vec = result[1]                # Optimized vector
optimal_metric = result[2]             # Optimal metric value

optimal_params = vector_to_param(optimal_vec, keys_list)
optimal_params = add_parameters(optimal_params)


println("-----------------------------------------------------")
println("-----------------------------------------------------")
println("------------------SIMULATION RESULTS-----------------")
println("-----------------------------------------------------")


cost(optimal_vec)

println("Optimal Parameters: $optimal_params")
println("Optimal Metric: $optimal_metric")
#println("Cg DIELECTRIC THICKNESS: ", optimal_params[:CgDielectricThichness])

circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, optimal_params, fixed_params)

S21, _, S11, _, S21phase = simulate_low_pump_power(sim_vars, circuit_temp, circuitdefs_temp)

maxS11 = maxS11val_BandFreq_FixFlux(S11, optimal_params, sim_vars)
println("Maximum S11: ", maxS11)


p1,p2,p3,p4 = plot_low_pump_power(S21, S11, S21phase, optimal_params, sim_vars)

p=plot(p1,p4,layout=(2,1), size=(600, 700))
display(p)


p, _, _ = plot_derivative_low_pump(S21phase, sim_vars, optimal_params)
display(p)

p_temp = simulate_and_plot(optimal_params, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)

println("Report completed")
display(p_temp)