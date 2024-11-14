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
    :smallJunctionArea => collect(2:0.5:5),                        
    :alphaSNAIL => collect(0.1:0.05:0.25),                         
    :LloadingCell => [1, 2, 3],                           
    :CgloadingCell => [1, 2],                             
    :criticalCurrentDensity => [0.2, 0.5],                
    :CgDielectricThichness => [10, 20,40,70]                    
)

nMacrocells_array = [3,20,40]
nMacrocells=20
nRandomPoints = 40

random_points=[]

for _ in 1:nRandomPoints

    params_temp = Dict(key => rand(values) for (key, values) in sim_params)
    
    #Adding important parameters
    params_temp[:N] = nMacrocells*params_temp[:loadingpitch] 
    params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
    params_temp[:CgAreaUNLoaded] = 150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])

    push!(random_points, params_temp)

end

for params_temp in random_points

    circuit_temp, circuitdefs_temp = create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params)

    start_time = now()  
    
    p1_temp, p2_temp, p3_temp, p4_temp = simulate_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp)
        
    p1p_temp, p2p_temp, p5_temp = simulate_at_fixed_flux(sim_vars, circuit_temp, circuitdefs_temp)

    p_temp = final_report(params_temp, sim_vars, fixed_params, p1_temp, p2_temp, p3_temp, p4_temp, p1p_temp, p2p_temp, p5_temp)

    end_time = now()

    display(p4_temp)

    println("Time taken: ", end_time - start_time)

end    





"""

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
        p1_temp, p2_temp, p3_temp, p4_temp = simulate_low_pump_power(params_temp, sim_vars, circuit_temp, circuitdefs_temp)
        
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

