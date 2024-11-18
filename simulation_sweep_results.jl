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
#using GaussianProcesses
using Distributions
using Surrogates

include("snail_circuit.jl")
include("simulation_black_box.jl")



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
    :fp =>  8.001 * 1e9,
    :Ip => 0.00001e-6,
    :IpGain => 0.8e-6,
    :IpSweep => (0.8:0.1:0.8) * 1e-6,
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
sim_params = Dict(
    :loadingpitch => [2, 3, 4],
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

function vector_to_param(vec, keys)

    keys_array = collect(keys)  # Convert KeySet to an array
    Dict(keys_array[i] => vec[i] for i in 1:length(keys_array))

end


function add_parameters(params_temp)

    #Adding important parameters
    params_temp[:N] = params_temp[:nMacrocells]*params_temp[:loadingpitch] 
    params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
    params_temp[:CgAreaUNLoaded] = 150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])

    return params_temp
end





# Define the objective function
function objective(vec) #vector passing
    
    # to add: funzione che da array passa a dictionary
    println("Vector:", vec)

    params_temp = vector_to_param(vec, keys(sim_params))
    println("Dictionary:", params_temp)

    params_temp = add_parameters(params_temp)
    println("Dictionary update:", params_temp)

    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)
    println("circuit created")
    
    alpha_wphalf, alpha_wp, _  = calculation_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp) 
    println("Metric calculated")

    return abs(alpha_wphalf - alpha_wp)

end



# Function to generate a single initial point as a vector of vectors
function generate_initial_point(params)
    return [rand(v) for v in values(params)]  # Correctly sampling from the values of the dictionary
end

# Function to generate n initial points
function generate_n_initial_points(n, params)
    return [generate_initial_point(params) for _ in 1:n]  # Generate n points
end


bounds = [extrema(v) for v in values(sim_params)]

n_initial_points = 4
initial_points = generate_n_initial_points(n_initial_points, sim_params)
println("initial_points: ", initial_points)

for p in initial_points
    println(p)
end 



initial_values = [objective(p) for p in initial_points]
println("initial_values: ", initial_values)

kernel = SquaredExponential()

# Perform Gaussian Process optimization
result = Surrogates.optimize(
    objective, 
    bounds, 
    Surrogates.GaussianProcess(kernel=kernel, noise_var=1.0), 
    maxiters = 100,
    initial_points = initial_points,
    initial_values = initial_values
)

keys_list = collect(keys(sim_params))  # Extract parameter names
optimal_vec = result[1]                # Optimized vector
optimal_metric = result[2]             # Optimal metric value

optimal_params = vector_to_param(optimal_vec, keys_list)

println("Optimal Parameters: $optimal_params")
println("Optimal Metric: $optimal_metric")











"""
function generate_random_point()
    
    return Dict(key => rand(values) for (key, values) in sim_params)

end



#create an array of random points
nRandomPoints = 3
random_points = []

for _ in 1:nRandomPoints

    params_temp=generate_random_point()

    push!(random_points, params_temp)

end

println(random_points[1])
"""


"""
wp_alpha_vec=[]
wphalf_alpha_vec=[]
lin_alpha_vec=[]

for params_temp in random_points

    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)

    start_time = now()  
    
    alpha_wphalf, alpha_wp, alpha_lin, p4_temp = simulate_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp)  
    
    push!(wphalf_alpha_vec, alpha_wphalf)
    push!(wp_alpha_vec, alpha_wp)
    push!(lin_alpha_vec, alpha_lin)


    p1p_temp, p2p_temp, p5_temp = simulate_at_fixed_flux(sim_vars, circuit_temp, circuitdefs_temp)
    p_temp = final_report(params_temp, sim_vars, fixed_params, p1_temp, p2_temp, p3_temp, p4_temp, p1p_temp, p2p_temp, p5_temp)


    end_time = now()

    display(p4_temp)

    println("Time taken: ", (end_time - start_time)/Second(1))

end    

delta_alpha_wp_vec =  abs.(wp_alpha_vec .- wphalf_alpha_vec)
println(delta_alpha_wp_vec)
delta_alpha_lin_vec =  abs.(lin_alpha_vec .- wphalf_alpha_vec)
println(delta_alpha_lin_vec)





function simulation_for_optmz(params_temp, JJSmallStd, JJBigStd, fixed_params, sim_vars)
    
    #Adding important parameters
    params_temp[:N] = nMacrocells*params_temp[:loadingpitch] 
    params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
    params_temp[:CgAreaUNLoaded] = 150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])

    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)
    
    alpha_wphalf, alpha_wp, _, _ = simulate_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp) 
    
    return alpha_wphalf, alpha_wp

end



function cost_func(params, JJSmallStd, JJBigStd, fixed_params, sim_vars) #, delta_alpha_lin_vec, threshold)

    alpha_wphalf, alpha_wp = simulation_for_optmz(params, JJSmallStd, JJBigStd, fixed_params, sim_vars)

    return abs(alpha_wphalf - alpha_wp)

end


initial_params = generate_random_point()

function sequential_minimization(params, nRandomPoints)

    for i in 1:nRandomPoints

        # Perform the optimization for the i-th element
        result = optimize(p -> cost_func(p, JJSmallStd, JJBigStd, fixed_params, sim_vars), params)

        # Update the parameters after each optimization
        params = result.minimizer
        
        # Display the progress and result for each element
        println("Optimal params for element dollaro i: ", params)
        println("Minimum cost for element dollaro i: ", result.minimum)
    end
    
    return params

end












#Valuation of the time as a function of the cells number

nMacrocells_array = [3,20,40,80,100,130,160,200,250]
loadingpitch = 3
time_array=[]

params_temp = Dict(key => rand(values) for (key, values) in sim_params)
params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
params_temp[:CgAreaUNLoaded] = 150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])

for nMacrocells in nMacrocells_array

    start_time = now() 

    params_temp[:N] = nMacrocells*loadingpitch
    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)
    p_temp = simulate_and_plot(params_temp, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)

    end_time = now() 

    push!(time_array, (end_time - start_time)/Second(1))

end

println(nMacrocells_array)
println(time_array)

plt = plot(nMacrocells_array, time_array,
    title="Time as a function of the number of cells", xlabel="nMacrocells", ylabel="time [s]")

display(plt)




--------------------------------------------------------------------------------------------------




# WARNING: for now choose only ONE parameter to iterate that YOU COMMENT below. In order to iterate it and fix all the others.

# Fix the variables that you don't want to iterate
fixed_vars = Dict(
    
    :loadingpitch => 2,
    :nMacrocells => 20,                #Number of macrocell composed by loadingpitch. The total number of cells N=nMacrocells*loadingpitch.
    :smallJunctionArea => 1,           #Tra 2 e 5 (meglio basso)   
    #:alphaSNAIL => 0.2,                #Tra 0.1 e 0.25            
    :LloadingCell => 2,                #Tra 1 e 6
    :CgloadingCell => 1,               #Tra 1 e 6
    :criticalCurrentDensity => 1,       #Tra 0.2 e 1         
    :CgDielectricThichness => 80         #Tra 10 e 70
)

# Combine fixed variables with the current value for simulation
params = merge(sim_params, fixed_vars)





#-----------------------------------------------------------------------------------
#-----------------------------STARTING SIMULATION-----------------------------------


for (param_name, param_values) in params
    
    if haskey(fixed_vars, param_name)
        continue
    end

    println("Simulating: ", param_name)  
    println("-----------------------------------------")

    # Loop through the values of the selected variable
    for value in param_values


        params_temp = merge(params, Dict(param_name => value))
        
        #Adding important parameters
        params_temp[:N] = params_temp[:nMacrocells]*params_temp[:loadingpitch] 
        params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
        params_temp[:CgAreaUNLoaded] = 150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])

        

        # Circuit creation
        circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)


        # Run the simulation and generate the plot

        # Unpack results with renamed variables
        p1_temp, p2_temp, p3_temp, p4_temp = plot_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp)
        
        p1p_temp, p2p_temp, p5_temp = simulate_at_fixed_flux(sim_vars, circuit_temp, circuitdefs_temp)

        p_temp = final_report(params_temp, sim_vars, fixed_params, p1_temp, p2_temp, p3_temp, p4_temp, p1p_temp, p2p_temp, p5_temp)

        #p_temp = simulate_and_plot(params_temp, sim_vars, fixed_params, circuit_temp, circuitdefs_temp)

    
        # Display the plot
        #display(p_temp)
        display(p4_temp)


     
        println("-----------------------------------------")

    end

end

"""