function vector_to_param(vec, keys)

    keys_array = collect(keys)  # Convert KeySet to an array
    Dict(keys_array[i] => vec[i] for i in 1:length(keys_array))

end


function add_parameters(params_temp)

    #Adding important parameters
    params_temp[:N] = params_temp[:nMacrocells]*params_temp[:loadingpitch] 
    params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
    params_temp[:CgAreaUNLoaded] = 150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])
    params_temp[:phidc] = find_flux_from_alpha(params_temp) 

    return params_temp
end


function find_n_initial_points(sim_params_space::Dict)

    parameter_lengths = [length(values) for values in values(sim_params_space)]
    n_initial_points = prod(parameter_lengths)
    
    return n_initial_points
end


# Function to generate a single initial point as a vector of vectors
function generate_initial_point(params)
    return [rand(v) for v in values(params)]  # Correctly sampling from the values of the dictionary
end


# Function to generate n initial points as a tuple with all float values
function generate_n_initial_random_points(n, params)
    # Generate `n` initial points as tuples with Float64 elements
    points = [Tuple(Float64.(generate_initial_point(params))) for _ in 1:n]
    return points
end


function generate_all_initial_points(params_space)
    # Extract the parameter value lists from the dictionary
    value_lists = values(params_space)
    
    # Initialize an empty set to store the generated points (sets do not allow duplicates)
    points_set = Set{Tuple{Float64, Vararg{Float64}}}()
    
    # Generate all combinations using nested for loops (generalized)
    for values in Iterators.product(value_lists...)
        # Convert each combination into a tuple of floats and add it to the set
        push!(points_set, tuple(map(float, values)...))
    end

    # Convert the set back to a list (array)
    return collect(points_set)
end

function save_points_to_file(points, filename)
    # Save the points as a vector to the file
    @save filename points
end

function load_points_from_file(filename)
    # Load the points vector from the file
    data = load(filename)
    return data["points"]  # Assuming the variable name is "points"
end


function simulation_time_estimation(n_initial_points, n_maxiters, n_num_new_samples)
    # Define the time per point in seconds
    time_per_point = 5.0 # seconds
    
    # Calculate total time for the simulation
    total_points = n_initial_points + n_maxiters * n_num_new_samples
    time_estimated = total_points * time_per_point  # in seconds
    
    # Breakdown time into days, hours, minutes, and seconds
    total_seconds = round(Int, time_estimated)
    days = div(total_seconds, 86400)
    hours = div(total_seconds % 86400, 3600)
    minutes = div(total_seconds % 3600, 60)
    seconds = total_seconds % 60
    
    # Calculate the finishing time
    current_time = Dates.now()  # Current date and time
    finish_time = current_time + Dates.Second(total_seconds)
    
    # Create formatted time strings
    formatted_estimation = "$(days)d $(hours)h $(minutes)m $(seconds)s"
    formatted_finish_time = string(finish_time) 
    
    return formatted_estimation, formatted_finish_time
end


function simulation_time(start_time)

    total_time = time() - start_time

    total_seconds = round(Int, total_time)
    days = div(total_seconds, 86400)
    hours = div(total_seconds % 86400, 3600)
    minutes = div(total_seconds % 3600, 60)
    seconds = total_seconds % 60

    formatted_time = "$(days)d $(hours)h $(minutes)m $(seconds)s"
    
    return formatted_time

end