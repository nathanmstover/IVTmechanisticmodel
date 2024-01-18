function getsupersaturation(parameters,Mgtotal,PPitotal,NTPtotal=1e-18,Buffertotal=1e-18)
    #Prevent divide by zero errors
    for species in [Mgtotal,PPitotal,NTPtotal,Buffertotal]
        if species<1e-18
            species = 1e-18
        end
    end
    #Calculate free species concentrations
    free = getfreeconcentrations(parameters, NTPtotal, Mgtotal, Buffertotal, PPitotal, 1e-18)
    #Calculate (ficticious) concentration of Magnesium Pyrophosphate
    (ionicspecies,_) = speciationmodel(parameters, NTPtotal, Mgtotal, Buffertotal, PPitotal, 1e-8, free...)
    Mg = ionicspecies[1]
    MgPPi = ionicspecies[3]
    Mg2PPi = ionicspecies[4]
    #Supersaturation
    return [Mg2PPi/p.Mg2PPi_eq, Mg2PPi, MgPPi, Mg]
end

function supersaturationprediction(model, parametersls, covariancematrix, inputs, n_mc_samples=10000, alpha = 0.05)
    n_outputs = 4
    mcensemble = zeros(n_outputs,n_mc_samples)
    mean = parametersls
    meanparams = fullparameterset(model,mean)
    d = MvNormal(mean, Hermitian(covariancematrix))
    for i in 1:n_mc_samples
        x = rand(d, 1)
        sampleparams = fullparameterset(model,x)
        sol = getsupersaturation(sampleparams,inputs...)
        mcensemble[:,i] = sol
    end

    maximumliklihood = getsupersaturation(meanparams,inputs...)
    lowerCB = [percentile(mcensemble[i,:],100*alpha/2) for i in 1:n_outputs]
    upperCB = [percentile(mcensemble[i,:],100*(1-alpha/2)) for i in 1:n_outputs]
    return maximumliklihood, lowerCB, upperCB
end