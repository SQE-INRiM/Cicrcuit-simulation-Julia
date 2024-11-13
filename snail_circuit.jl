#-------------------------------------CREATION OF SNAIL CIRCUIT-------------------------------------------

function create_circuit(JJSmallStd, JJBigStd, params_temp, fixed_params, verbose = false)
    
    #Circuit definitions----------------------------------------------------

    @variables Rleft Rright Cg Lj Cj alpha Lloading Cgloading M Ladd Lf Cf Lg
    

    circuitdefs = Dict(
        alpha => params_temp[:alphaSNAIL],                       # Directly access params_temp
        Lloading => params_temp[:LloadingCell],                  # Directly access fixed_params
        Cgloading => params_temp[:CgloadingCell],                # Directly access fixed_params
        
        Lj => IctoLj(params_temp[:smallJunctionArea] * params_temp[:criticalCurrentDensity] * 1e-6),  # H
        Cg => params_temp[:CgAreaUNLoaded] * params_temp[:CgDensity] / (1 + im * fixed_params[:tandeltaCgDielectric]),  # F
        Cj => params_temp[:smallJunctionArea] * fixed_params[:JosephsonCapacitanceDensity] * 1e-15,  # F

        #circut parameters
        Rleft => 50.0,            # Ohm
        Rright => 50.0,           # Ohm
        Ladd => 70.0e-15,         # Henry - loop inductance
        Lg => 20.0e-9,            # Henry - geometrical inductance
        Lf => 190.0e-12,          # Henry - Flux line inductors
        Cf => 0.076e-12,          # Farad - Flux line capacitors
        M => 0.999                # inverse inductance matrix with K<1.0
    )



    #Circuit creation-------------------------------------------------------------

    circuit = Tuple{String,String,String,Num}[]

    # Port on the input side of the AC line
    push!(circuit,("P$(1)_$(0)","1","0",1))
    push!(circuit,("R$(1)_$(0)","1","0",Rleft))

    if verbose
        println("Port 1 created")
        println("AC Node $(0) created")
        println("AC Node $(1) created")
    end

    rngSmall1 = MersenneTwister(1);
    randomSeedSmall1 = 1 + JJSmallStd*randn(rngSmall1, Float64)
    rngBig1 = MersenneTwister(1);
    randomSeedBig1 = 1 + JJBigStd*randn(rngBig1, Float64)

    #AC line---------------------------------------------------
    #first cell------------------------------------------------

    #first half capacitance to ground
    push!(circuit,("C$(0)_$(1)","0","1",Cg/2))
    push!(circuit,("Lj$(1)_$(2)","1","2",alpha*Lj*randomSeedBig1))                          # First big JJ inductance 
    push!(circuit,("C$(1)_$(2)","1","2",Cj/(alpha*randomSeedBig1)))                         # First big JJ capacitance 
    push!(circuit,("Lj$(2)_$(3)","2","3",alpha*Lj*randomSeedBig1))                          # Second big JJ inductance 
    push!(circuit,("C$(2)_$(3)","2","3",Cj/(alpha*randomSeedBig1)))                         # Second big JJ capacitance 
    push!(circuit,("Lj$(3)_$(5)","3","5",alpha*Lj*randomSeedBig1))                          # Third big JJ inductance 
    push!(circuit,("C$(3)_$(5)","3","5",Cj/(alpha*randomSeedBig1)))                         # Third big JJ capacitance 

    push!(circuit,("Lj$(1)_$(4)","1","4",Lj*randomSeedSmall1))                              # Small JJ inductance 
    push!(circuit,("C$(1)_$(4)","1","4",Cj/(randomSeedSmall1)))                       # Small JJ capacitance 

    push!(circuit,("L$(4)_$(5)","4","5",Ladd))                                              # Loop inductance for flux bias

    if verbose
        println("1: UNLoaded cell")
        println("AC Node 2 created")
        println("AC Node 3 created")
    end

    nodePerCell=fixed_params[:nodePerCell]
    N=params_temp[:N]
    loadingpitch=params_temp[:loadingpitch]

    j=nodePerCell+1

    for i = 2:N
        
        local rngSmall = MersenneTwister(i+1);
        local randomSeedSmall1 = 1+JJSmallStd*randn(rngSmall, Float64)
        local rngBig1 = MersenneTwister((i+1)*j+1);
        local randomSeedBig1 = 1+JJBigStd*randn(rngBig1, Float64)


        if mod(i,loadingpitch)+1 == loadingpitch÷2            
            #make the loaded cell
            push!(circuit,("C$(0)_$(j)","$(0)","$(j)",Cgloading*Cg))                                                    #capacitance to ground
            push!(circuit,("Lj$(j+0)_$(j+1)","$(j+0)","$(j+1)",Lloading*alpha*Lj*randomSeedBig1))                           # First big JJ inductance 
            push!(circuit,("C$(j+0)_$(j+1)","$(j+0)","$(j+1)",((1/Lloading)*Cj)/(alpha*randomSeedBig1)))                    # First big JJ capacitance 
            push!(circuit,("Lj$(j+1)_$(j+2)","$(j+1)","$(j+2)",Lloading*alpha*Lj*randomSeedBig1))                       # Second big JJ inductance 
            push!(circuit,("C$(j+1)_$(j+2)","$(j+1)","$(j+2)",((1/Lloading)*Cj)/(alpha*randomSeedBig1)))                # Second big JJ capacitance 
            push!(circuit,("Lj$(j+2)_$(j+4)","$(j+2)","$(j+4)",Lloading*alpha*Lj*randomSeedBig1))                       # Third big JJ inductance 
            push!(circuit,("C$(j+2)_$(j+4)","$(j+2)","$(j+4)",((1/Lloading)*Cj)/(alpha*randomSeedBig1)))                # Third big JJ capacitance 

            push!(circuit,("Lj$(j)_$(j+3)","$(j)","$(j+3)",Lloading*Lj*randomSeedSmall1))                               # Small JJ inductance 
            push!(circuit,("C$(j)_$(j+3)","$(j)","$(j+3)",((1/Lloading)*Cj)/randomSeedSmall1))                          # Small JJ capacitance 
        
            push!(circuit,("L$(j+3)_$(j+4)","$(j+3)","$(j+4)",Ladd))     # Loop inductance for flux bias
            
            if verbose
                println("$(i): Loaded cell")                                                # Loop inductance for flux bias
            end

        else
            # make the unloaded cell
            push!(circuit,("C$(0)_$(j)","$(0)","$(j)",Cg))                                                              #capacitance to ground
            push!(circuit,("Lj$(j+0)_$(j+1)","$(j+0)","$(j+1)",alpha*Lj*randomSeedBig1))                                    # First big JJ inductance 
            push!(circuit,("C$(j+0)_$(j+1)","$(j+0)","$(j+1)",Cj/(alpha*randomSeedBig1)))                                   # First big JJ capacitance 
            push!(circuit,("Lj$(j+1)_$(j+2)","$(j+1)","$(j+2)",alpha*Lj*randomSeedBig1))                                # Second big JJ inductance 
            push!(circuit,("C$(j+1)_$(j+2)","$(j+1)","$(j+2)",Cj/(alpha*randomSeedBig1)))                               # Second big JJ capacitance 
            push!(circuit,("Lj$(j+2)_$(j+4)","$(j+2)","$(j+4)",alpha*Lj*randomSeedBig1))                                # Third big JJ inductance 
            push!(circuit,("C$(j+2)_$(j+4)","$(j+2)","$(j+4)",Cj/(alpha*randomSeedBig1)))                               # Third big JJ capacitance 

            push!(circuit,("Lj$(j)_$(j+3)","$(j)","$(j+3)",Lj*randomSeedSmall1))                                        # Small JJ inductance 
            push!(circuit,("C$(j)_$(j+3)","$(j)","$(j+3)",Cj/randomSeedSmall1))                                         # Small JJ capacitance 
        
            push!(circuit,("L$(j+3)_$(j+4)","$(j+3)","$(j+4)",Ladd))       # Loop inductance for flux bias
            
            if verbose
                println("$(i): UNLoaded cell")                                                # Loop inductance for flux bias
            end

        end

        if verbose
            println("AC Node j+1=$(j+1) created")
            println("AC Node j+2=$(j+2) created")
        end
        
        # increment the index
        j = j+nodePerCell

    end

    #last cell
    push!(circuit,("C$(0)_$(j)","$(0)","$(j)",Cg/2))
    push!(circuit,("R$(0)_$(j)","$(0)","$(j)",Rright))

    #AC port on the output side
    push!(circuit,("P$(0)_$(j)","$(0)","$(j)",2))
    
    if verbose
        println("Port 2 created")
    end
    
    #END AC line--------------------------------------------------------------------------------------


    #DC line------------------------------------------------------------------------------------------

    # port on the input side of the DC line
    dcOffs = nodePerCell*N+2+1
    push!(circuit,("P$(dcOffs+1)_$(0)","$(dcOffs+1)","$(0)",3))
    push!(circuit,("R$(dcOffs+1)_$(0)","$(dcOffs+1)","$(0)",Rleft))

    if verbose
        println("Port 3 created")
        println("DC Node $(dcOffs) created")
        println("DC Node $(dcOffs+1) created")
    end

    
    #first cell---------------------------------------------------------------------------------------
    push!(circuit,("C$(dcOffs+1)_$(0)","$(dcOffs+1)","$(0)",Cf/2))                                          #DC line capacitance in the first cell
    push!(circuit,("L$(dcOffs+1)_$(dcOffs+2)","$(dcOffs+1)","$(dcOffs+2)",Lf))                              #DC line inductance in the first cell
    push!(circuit,("K$(1)_$(1)","L$(4)_$(5)","L$(dcOffs+1)_$(dcOffs+2)",M))                          #mutual inductance between loop inductance and DC line
    
    if verbose
        println("L$(4)_$(5) mutually connected to L$(dcOffs+1)_$(dcOffs+2)")
        println("DC Node dcOffs+1+1=$(dcOffs+1+1) created")
        println("K$(1)_$(1) mutual coupling created")
    end
    
    for i = 2:N
        #DC line
        push!(circuit,("C$(dcOffs+i)_$(0)","$(dcOffs+i)","$(0)",Cf))                                        #DC line capacitance
        push!(circuit,("L$(dcOffs+i)_$(dcOffs+i+1)","$(dcOffs+i)","$(dcOffs+i+1)",Lf))                      #DC line inductance
        #AC-DC mutual coupling
        push!(circuit,("K$(i)_$(i)","L$(i*nodePerCell)_$(i*nodePerCell+1)","L$(dcOffs+i)_$(dcOffs+i+1)",M)) #mutual inductance between loop inductance and DC line (equal for each cell)
        
        if verbose
            println("$(2*(i % 2-0.5)): mutual inductance sign")
            println("L$(i*nodePerCell)_$(i*nodePerCell+1) mutually connected to L$(dcOffs+i)_$(dcOffs+i+1)")
            println("DC Node dcOffs+i+1=$(dcOffs+i+1) created")
            println("K$(i)_$(i) mutual coupling created")
        end
    
    end

    #DC port on the output side
    push!(circuit,("P$(dcOffs+1+N)_$(0)","$(dcOffs+1+N)","$(0)",4))
    push!(circuit,("L$(dcOffs+1+N)_$(0)","$(dcOffs+1+N)","$(0)",Lg))
    push!(circuit,("C$(dcOffs+1+N)_$(0)","$(dcOffs+1+N)","$(0)",Cf/2))    
    push!(circuit,("R$(dcOffs+1+N)_$(0)","$(dcOffs+1+N)","$(0)",Rright))

    if verbose
        println("Port 4 created")
        println("$(dcOffs+1+N)")
    end

    #END DC line--------------------------------------------------------------------------------------


    return circuit, circuitdefs

end


