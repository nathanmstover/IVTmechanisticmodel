#Global variables for species indexing
RNAindex = 2
PPiindex = 3
Mgindex = 8
#Global variable defines 95% prediction intervals (99% -> alpha = 0.01)
alpha = 0.05

"""
    autodiffgradient(speciesindex, params, inputs; kargs...)

Use forward mode autodifferentiation to genrate gradient of output (speciesindex) with respect to parameters.
"""
function autodiffgradient(speciesindex, params, inputs; kargs...)
    func = x -> (runDAE_batch(x, inputs; kargs...).u[end])[speciesindex]
    autodiffgrad = ForwardDiff.gradient(func, params)
    funceval = func(params)
    return funceval,autodiffgrad
end

"""
    PPititrationsensitivity(model,data,fitparameters,log=true)

Generate sensitivity matrix of Akama Mg2PPi solubility data. Either outputs derivative of output with respect to paramters, or derivative of the log of the output with respect to parameters.
"""
function PPititrationsensitivity(model,data,fitparameters,log=true)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data.Mgconcentrationinputs)[1]
    numberofparameters = length(iterationparameters)
    numberoftimepoints = length(data.timeinputs[1,:])
    
    sensitivitymatrix = zeros(numberofparameters,numberofindependentconditions)
    
    for i in 1:numberofindependentconditions
            
            NTPconc = data.Mgconcentrationinputs[i,1]
            PPiconc = data.Mgconcentrationinputs[i,2] 

            inputs =(T7RNAP = 0, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = 0.004, Buffer = 0.040, DNA = 1e-10, final_time = 24)
      
            eval,autodiffgrad = autodiffgradient(Mgindex,iterationparameters,inputs; saveevery=false, PPi = PPiconc)

            if log 
                evaluationvalue = eval
            else
                evaluationvalue = 1
            end
            sensitivitymatrix[:,i] = Base.log(10)*iterationparameters .* autodiffgrad ./ evaluationvalue
    end
    return sensitivitymatrix
end

"""
    initialratesensitivity(model,data,fitparameters,log=true)

Generate sensitivity matrix of initial Akama transcription rate data. Either outputs derivative of output with respect to paramters, or derivative of the log of the output with respect to parameters.
"""
function initialratesensitivity(model,data,fitparameters,log=true)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data.initialrateconcentrationinputs)[1]
    numberofparameters = length(iterationparameters)
    
    sensitivitymatrix = zeros(numberofparameters,numberofindependentconditions)
    
    for i in 1:numberofindependentconditions
            NTPconc = data.initialrateconcentrationinputs[i,1]
            Mgconc = data.initialrateconcentrationinputs[i,2]

            inputs =(T7RNAP = 1e-7, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = Mgconc, Buffer = 0.040, DNA = 7.4, final_time = 0.08333333333)
            eval,autodiffgrad = autodiffgradient(RNAindex,iterationparameters,inputs; saveevery=false)

            if log 
                evaluationvalue = eval
            else
                evaluationvalue = 1
            end
            sensitivitymatrix[:,i]= Base.log(10)*iterationparameters .* autodiffgrad ./ evaluationvalue
    end
    return sensitivitymatrix
end

"""
    sensitivitymatricies(model,data,fitparameters,log=true)

Generate sensitivity matrix of solution trajectories Akama data, for both RNA and PPi trajectories. Either outputs derivative of output with respect to paramters, or derivative of the log of the output with respect to parameters.
"""
function sensitivitymatricies(model,data,fitparameters,log=true)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data.concentrationinputs)[1]
    numberofparameters = length(iterationparameters)
    numberoftimepoints = length(data.timeinputs[1,:])
    RNAsensitivitymatrix = zeros(length(iterationparameters),numberofindependentconditions*numberoftimepoints)
    PPisensitivitymatrix = zeros(length(iterationparameters),numberofindependentconditions*numberoftimepoints)

    for i in 1:numberofindependentconditions
        for (j,time) in enumerate(data.timeinputs[1,:])
            
            T7 = data.concentrationinputs[i,1]
            NTPtot = data.concentrationinputs[i,2]
            Mg = data.concentrationinputs[i,3]
            
            inputs = (T7RNAP = T7, ATP = NTPtot/4,UTP = NTPtot/4,CTP = NTPtot/4,GTP = NTPtot/4, Mg = Mg, Buffer = 0.040, DNA = 7.4, final_time = time)
            
            RNAeval, RNAautodiffgrad = autodiffgradient(RNAindex,iterationparameters,inputs; saveevery=false)
            if log == false
                RNAeval = 1
            end
            RNAsensitivitymatrix[:,(i-1)*numberoftimepoints+j]= Base.log(10)*iterationparameters .* RNAautodiffgrad ./ RNAeval
            
            PPieval, PPiautodiffgrad = autodiffgradient(PPiindex,iterationparameters,inputs; saveevery=false)
            if log != true
                PPieval = 1
            end
            PPisensitivitymatrix[:,(i-1)*numberoftimepoints+j]= Base.log(10)*iterationparameters .* PPiautodiffgrad ./ PPieval  
        end
    end
    return RNAsensitivitymatrix,PPisensitivitymatrix
end

"""
    sensitivitymatricies(model,data,fitparameters,log=true)

Generate sensitivity matrix of solution trajectories Akama data, for both RNA and PPi trajectories. Either outputs derivative of output with respect to paramters, or derivative of the log of the output with respect to parameters.
"""
function customsensitivitymatrix(model,data,fitparameters,log = true)
    iterationparameters = fullparameterset(model,fitparameters)
    numberofindependentconditions = size(data)[1]
    sensitivitylist = []
    stdevs = []

    for i in 1:numberofindependentconditions
        T7RNAP = data[i,3]/1e9
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

        RNAconfinput = parseexceltuple(data[i,20])

        if speciesindex == 2
            scalingfactor = 1e6
        else
            scalingfactor = 1e3
        end

        #Can either input uncertainty as scalar or vector
        if length(RNAconfinput) == 1
            RNAconfint = RNAconfinput .* ones(length(times))
        else
            RNAconfint = RNAconfinput
        end

        stoich = SVector(N_A, N_U, N_C, N_G)

        for (j,time) in enumerate(times)

            inputs =(T7RNAP = T7RNAP, ATP = ATP,UTP = UTP,CTP = CTP,GTP = GTP, Mg = Mg, Buffer = Buffer, DNA = DNA, final_time = time)
            eval, autodiffgrad = autodiffgradient(speciesindex,iterationparameters,inputs; saveevery=false, Cap = Cap, PPiase = PPiase, stoich = stoich)
            
            if log == false
                RNAeval = 1
            end
            sensitivity = Base.log(10)*iterationparameters .* autodiffgrad ./ RNAeval
            append!(sensitivitylist, [sensitivity])
            append!(stdevs,RNAconfint[j]/scalingfactor)
        end
    end
    sensitivitymatrix = hcat(sensitivitylist...)
    return sensitivitymatrix, stdevs
end