RNAindex = 2
PPiindex = 3
Mgindex = 8



function printresidual(model, data, x; customfile = false, customfilename = "")
    print("Total Residual:                        ")
    if customfile
        println(string(round(customresidualeval(model,data,customfilename,x),digits = 1)))
    else
        println(string(round(akamaresidual(model,data,x),digits = 1)))
    end
    println("Components of residual:")
    println("Concentration Trajectories (Figure 2): "*string(round(trajectorydataresidual(model,data,x),digits = 1)))
    println("Initial Reaction Rate (Figure 3A):     "*string(round(initialrateresidual(model,data,x),digits = 1)))
    println("Mg2PPi solubility (Figure 3B):         "*string(round(mg2ppisolubilityresidual(model,data,x),digits = 1)))
    println("Parameter priors:                      "*string(round(parameterresidual(model,x),digits = 1)))
    if customfile
        println("Custom Data:                           "*string(round(customresidual(model, customfilename, x),digits = 1)))
    end
end

"""
    akamaresidual(model,data,x)

Take model and candidate parameter vector, return residual for only Akama dataset.
""" 
function akamaresidual(model,data,x)
    residual = 0
    residual+=mg2ppisolubilityresidual(model,data,x)
    residual+=trajectorydataresidual(model,data,x)
    residual+=parameterresidual(model,x)
    residual+=initialrateresidual(model,data,x)
    return residual
end

"""
    customresidualeval(model,data,filename,x)

Takes csv filename, total residual for data file plus Akama data.
""" 
function customresidualeval(model,data,filename,x)
    residual = 0
    residual+=akamaresidual(model,data,x)
    residual+=customresidual(model, filename, x)
    return residual
end

"""
    customresidual(model, filename, fitparameters)

Takes csv filename, return residual for data in file.
""" 
function customresidual(model, filename, fitparameters)
    df = CSV.read(filename, DataFrame)
    data = Matrix(df)
    params = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data)[1]
    loss = 0
    for i in 1:numberofindependentconditions
        times = parseexceltuple(data[i,2])
        T7RNAP = data[i,3]/1e9#(8.41/6)*
        DNA = data[i,4]
        ATP = data[i,5]/1000
        UTP = data[i,6]/1000
        CTP = data[i,7]/1000
        GTP = data[i,8]/1000
        Cap = data[i,9]/1000

        Mg = data[i,10]/1000
        PPiase = data[i,11]
        Buffer = data[i,12]/1000

        N_A = data[i,15]
        N_U = data[i,16]
        N_C = data[i,17]
        N_G = data[i,18]
        
        speciesindex = data[i,1]
        times = parseexceltuple(data[i,2])
        RNAyields = parseexceltuple(data[i,19])

        RNAconfinput = parseexceltuple(data[i,20])

        if speciesindex == 2
            scalingfactor = 1e6
        else
            scalingfactor = 1e3
        end

        #Can either input uncertainty as scalar or vector
        if length(RNAconfinput) == 1
            RNAconfint = RNAconfinput[1] .* ones(length(RNAyields))
        else
            RNAconfint = RNAconfinput
        end
        stoich = SVector(N_A, N_U, N_C, N_G)

        inputs =(T7RNAP = T7RNAP, ATP = ATP,UTP = UTP,CTP = CTP,GTP = GTP, Mg = Mg, Buffer = Buffer, DNA = DNA, final_time = maximum(times))
        sol = runDAE_batch(params, inputs, PPiase = PPiase, stoich = stoich, Cap = Cap)

        RNApredictions = scalingfactor*sol(times)[speciesindex,:]
        loss+=norm((RNAyields .- RNApredictions) ./RNAconfint)^2
    end
    return loss
end

"""
    initialrateresidual(model, data, fitparameters)

Take model and candidate parameter vector, return residual for Akama initial rate dataset.
""" 
function initialrateresidual(model, data, fitparameters)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data.initialrateconcentrationinputs)[1]
    numberofparameters = length(iterationparameters)
    residual = 0
    for i in 1:numberofindependentconditions
        NTPconc = data.initialrateconcentrationinputs[i,1]
        Mgconc = data.initialrateconcentrationinputs[i,2]
        RNAyield = data.initialrateyieldoutputs[i]
        RNAyieldstdev = data.initialrateyieldstdev[i]
        inputs =(T7RNAP = 1e-7, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = Mgconc, Buffer = 0.040, DNA = 7.4, final_time = 0.08333333333)
        sol = runDAE_batch(iterationparameters, inputs; saveevery=false)
        modelRNAyield = sol.u[end][RNAindex]
        r=(modelRNAyield-RNAyield)^2/(RNAyieldstdev)^2
        residual += r
    end
    return residual
end

"""
    trajectorydataresidual(model, data, fitparameters)

Take model and candidate parameter vector, return residual for Akama concentration trajectory dataset.
""" 
function trajectorydataresidual(model, data, fitparameters)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data.concentrationinputs)[1]
    numberofparameters = length(iterationparameters)
    residual = 0
    for i in 1:numberofindependentconditions
        RNAPconc = data.concentrationinputs[i,1]
        NTPconc = data.concentrationinputs[i,2]
        Mgconc = data.concentrationinputs[i,3]

        inputs =(T7RNAP = RNAPconc, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = Mgconc, Buffer = 0.040, DNA = 7.4, final_time = 2)
        sol = runDAE_batch(iterationparameters, inputs)

        timepointpredictions = sol(data.timeinputs[i,:]).u
        timepointRNA = [i[RNAindex] for i in timepointpredictions]
        timepointPPi = [i[PPiindex] for i in timepointpredictions]
        r=norm((timepointRNA .- data.RNAyieldoutputs[i,:]) ./data.RNAstdev[i,:])^2 + norm((timepointPPi .- data.PPiyieldoutputs[i,:]) ./data.PPistdev[i,:] )^2
        residual += r
    end
    return residual
end

"""
    mg2ppisolubilityresidual(model, data, fitparameters)

Take model and candidate parameter vector, return for Akama Mg2PPi solubility dataset.
""" 
function mg2ppisolubilityresidual(model, data, fitparameters)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data.Mgconcentrationinputs)[1]
    numberofparameters = length(iterationparameters)
    residual = 0
    for i in 1:numberofindependentconditions


        NTPconc = data.Mgconcentrationinputs[i,1]
        PPiconc = data.Mgconcentrationinputs[i,2]
        Mgconc = data.Mgoutputs[i]
        Mgconcstdev = data.Mgstdev[i]

        inputs = (T7RNAP = 0, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = 0.004, Buffer = 0.040, DNA = 1e-8, final_time = 24)
        sol = runDAE_batch(iterationparameters, inputs; saveevery=false, PPi = PPiconc);
        modelMgoutput = sol.u[end][Mgindex]
        r= (modelMgoutput-Mgconc)^2/(Mgconcstdev)^2
        residual += r

    end
    return residual
end

"""
    parameterresidual(model,fitparameters)

Take model and candidate parameter vector, return residual representing deviation from Bayesian prior.
""" 
function parameterresidual(model,fitparameters)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofparameters = length(iterationparameters)
    residual = 0    
    for j in 1:numberofparameters
        if model.parameters[j].hasprior
            residual+= ((log10(model.parameters[j].value)-log10(iterationparameters[j]))^2)/(model.parameters[j].variance)^2
        end
    end
    return residual
end

