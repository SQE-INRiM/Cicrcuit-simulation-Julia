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
    :JosephsonCapacitanceDensity => 45, # Unit value
    :CgDielectricK => 9.6,             # Dielectric constant
    :tandeltaCgDielectric => 2.1e-3,   # Loss tangent of dielectric
)

#SIMULATION PARAMETERS

sim_vars = Dict(
    
    :fs => (0.1:0.1:30.0) * 1e9,
    :fp =>  14.001 * 1e9,
    :Ip => 0.00001e-6,
    :IpGain => 0.8e-6,
    :IpSweep => (0.8:0.1:0.8) * 1e-6, #0.8 iniziale
    :phidc => 0.38,
    :phidcSweep => (0:0.01:0.5),
    :Npumpharmonics => (8,),
    :Nmodulationharmonics => (4,),
    :Niterations => 500
)

sim_vars[:ws] = 2 * pi * sim_vars[:fs]
sim_vars[:wp] = (2 * pi * sim_vars[:fp],)


#--------------------------------------------------------------------------------


#PARAMETER SCAN

JJSmallStd = 0.05               #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread
JJBigStd = 0.05                 #Tra 0.05 e 0.2             # =0 -> perfect fab , 0.1 -> 10% spread

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
#-----------------------------STARTING SIMULATION-----------------------------------
#-----------------------------------------------------------------------------------


# Define the cost function
function cost(vec) #vector passing
    
    println("-----------------------------------------------------")
    println("Vector: ", vec)

    params_temp = vector_to_param(vec, keys(sim_params_space))
    #println("Dictionary: ", params_temp)

    params_temp = add_parameters(params_temp)
    println("Dictionary update: ", params_temp)

    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)
    println("Circuit created")

    S21, S12, S11, S22, S21phase = simulate_low_pump_power(sim_vars, circuit_temp, circuitdefs_temp)
    println("S parameters calculated")

    maxS11value = maxS11val_BandFreq_FixFlux(S11, params_temp, sim_vars)
    println("Maxim value of S11: ", maxS11value)

    alpha_wphalf, alpha_wp, alpha_lin  = calculation_low_pump_power(S21phase, params_temp, sim_vars)

    println("Metric calculated")

    delta_alpha_wp=abs(alpha_wp - alpha_lin)
    delta_alpha_wphalf=abs(alpha_wphalf - alpha_lin)


    metric=abs(delta_alpha_wp-delta_alpha_wphalf) * (1/((abs(delta_alpha_wp*delta_alpha_wphalf))^(1/2)))
    
    println("Value: ", metric)

    println("-----------------------------------------------------")

    if maxS11value<-20                      #dBm

        return metric

    else

        return exp(-metric)

    end

    #p_temp = simulate_and_plot(params_temp, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)
    #display(p_temp)

end




println("START")
#Return upper and lower boundaries
bounds = [extrema(v) for v in values(sim_params_space)]
println("Bounds: ", bounds)

lb = [b[1] for b in bounds]
ub = [b[2] for b in bounds]
println("Lower bounds:", lb)
println("Upper bounds:", ub)



n_initial_points = 10
initial_points = generate_n_initial_points(n_initial_points, sim_params_space)
println("Initial points: ", initial_points)

for p in initial_points
    println("Single point: ", p)
end 

initial_values = [cost(p) for p in initial_points]
println("Initial values: ", initial_values)


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
    maxiters = 10,           # Number of interactions. Incresing maxiters: Leads to a longer optimization process with potentially better solutions but at the cost of more time.
    num_new_samples = 10     # Number of point generated for every single interaction. Incresing num_new_samples: Allows each iteration to consider a broader range of candidate points, 
                            # improving the chance of finding a good solution early but also increasing the computational cost per iteration.
)


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





println("-----------------------------------------------------")
println("------------------SIMULATION RESULTS-----------------")
println("-----------------------------------------------------")
println("-----------------------------------------------------")


optimal_params = add_parameters(optimal_params)

println("Optimal Parameters: $optimal_params")
println("Optimal Metric: $optimal_metric")

circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, optimal_params, fixed_params)

p_temp = simulate_and_plot(optimal_params, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)

println("Report completed")
display(p_temp)