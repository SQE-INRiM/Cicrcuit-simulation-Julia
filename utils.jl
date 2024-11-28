function vector_to_param(vec, keys)

    keys_array = collect(keys)  # Convert KeySet to an array
    Dict(keys_array[i] => vec[i] for i in 1:length(keys_array))

end


function add_parameters(params_temp)

    #Adding important parameters
    params_temp[:N] = params_temp[:nMacrocells]*params_temp[:loadingpitch] 
    params_temp[:CgDensity] = (fixed_params[:CgDielectricK] * 8.854e-12) / (1e12 * params_temp[:CgDielectricThichness] * 1e-9)
    params_temp[:CgAreaUNLoaded] = 200 #150 + 20 * (params_temp[:smallJunctionArea] / params_temp[:alphaSNAIL])
    params_temp[:phidc] = find_flux_from_alpha(params_temp) 

    return params_temp
end


# Function to generate a single initial point as a vector of vectors
function generate_initial_point(params)
    return [rand(v) for v in values(params)]  # Correctly sampling from the values of the dictionary
end


# Function to generate n initial points as a tuple with all float values
function generate_n_initial_points(n, params)
    # Generate `n` initial points as tuples with Float64 elements
    points = [Tuple(Float64.(generate_initial_point(params))) for _ in 1:n]
    return points
end
