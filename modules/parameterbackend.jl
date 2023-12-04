"""
    IVTmodel

List of parameter objects."""
struct IVTmodel
    parameters
end

"""
    Parameter

Collects information of name, fitting status, prior (hasprior, prior value, and prior variance), and upper and lower bound of parameter for optimization."""
struct Parameter
    name
    isfitted
    hasprior
    value
    variance
    upperbound
    lowerbound
end

"""
    parametermask(model)

Take model, return bitvector representing if parameters are fitted or fixed.
"""
function parametermask(model)
    mask = zeros(length(model.parameters))
    for i in 1:length(model.parameters)
        if model.parameters[i].isfitted
            mask[i] = 1
        end
    end
    return mask      
end

Parameter(name,value) = Parameter(name,false,false,value,0,0,0)
Parameter(name,value,upperbound,lowerbound) = Parameter(name,true,false,value,0,upperbound,lowerbound)
Parameter(name,value,upperbound,lowerbound,standarddev) = Parameter(name,true,true,value,standarddev,upperbound,lowerbound)

"""
    fullparameterset(model,roundparameters)

Take model and list of fitted parameter values, return ComponentArray of all parameters for use in model evaluation.
"""
function fullparameterset(model,roundparameters)
    namesofvalues::Vector{String} = []
    parametervalues = []
    fittedparametercounter = 1
    for i in model.parameters
        append!(namesofvalues,[i.name])
        if i.isfitted
            append!(parametervalues,[10^(roundparameters[fittedparametercounter])])
            fittedparametercounter+=1
        else
            append!(parametervalues,[i.value])  
        end
    end
    nt = ComponentArray(namedtuple(namesofvalues, parametervalues))
    return nt
end

addtoinputlist!(list,model) = append!(list,[model])

"""
    setupmodel()

Generate model with parameter names, fitting status, and prior distrubution.
"""
function setupmodel()
    #This is our "interface": each line here represents adding a parameter and some associated information
    #format is Parameter(name, bayesian prior, lower bound (useful for optimization), upper bound, bayesian stdev)
    modelinputparameters = []
    #Parameters for transcription rate
    addtoinputlist!(modelinputparameters,Parameter("k_i",936,1,50000,0.3))
    addtoinputlist!(modelinputparameters,Parameter("k_e",5.3e5,1e4,15e5,0.15))
    addtoinputlist!(modelinputparameters,Parameter("k_off",4320,400,74000,0.25))
    addtoinputlist!(modelinputparameters,Parameter("k_on",204,20,2100,0.05))
    addtoinputlist!(modelinputparameters,Parameter("K_1",10^(-3.692),0.00001,0.1))
    addtoinputlist!(modelinputparameters,Parameter("K_2",10^(-3.80),0.00001,0.05))
    addtoinputlist!(modelinputparameters,Parameter("Ki_PPi",0.0002,0.00001,0.05,0.3))

    # #Parameters for Pyrophosphatase activity
    addtoinputlist!(modelinputparameters,Parameter("kPPiase",1,))
    addtoinputlist!(modelinputparameters,Parameter("KMPPiase",0.000214))
    # addtoinputlist!(modelinputparameters,Parameter("kPPiase",1,0.1,10,0.08))
    # addtoinputlist!(modelinputparameters,Parameter("KMPPiase",0.000214,0.0000214,0.00214,0.08))

    #Parameters for solid formation
    addtoinputlist!(modelinputparameters,Parameter("k_precip",10^(-0.26),1e-8,1e6))
    addtoinputlist!(modelinputparameters,Parameter("B",13.48,1e-5,90,0.1))
    addtoinputlist!(modelinputparameters,Parameter("k_d",10^(4.47),1e-2,1e12))

    #Parameters for equilibria
    addtoinputlist!(modelinputparameters,Parameter("K_HNTP",10^(6.91),10^(4.91),10^(8.91),0.02))#
    addtoinputlist!(modelinputparameters,Parameter("K_HMgNTP",10^(2.08),10^(0.08),10^(4.08),0.05))#
    addtoinputlist!(modelinputparameters,Parameter("K_HPPi",10^(9.02),10^(8.02),10^(10.01),0.05))
    addtoinputlist!(modelinputparameters,Parameter("K_HMgPPi",10^(3.32),10^(1.32),10^(4.32),0.05))      
    addtoinputlist!(modelinputparameters,Parameter("K_H2PPi",10^(6.26),10^(4.26),10^(8.26),0.1))#
    addtoinputlist!(modelinputparameters,Parameter("K_H2MgPPi",10^(2.11),10^(0.11),10^(4.11),0.1))#
    addtoinputlist!(modelinputparameters,Parameter("K_MgNTP",10^(4.54),10^(3.0),10^(6),0.58))
    addtoinputlist!(modelinputparameters,Parameter("K_Mg2NTP",10^(1.77),10^(0.01),10^(3.5),0.58))
    addtoinputlist!(modelinputparameters,Parameter("K_MgPPi",10^(4.80),10^(3.80),10^(7),0.58))
    addtoinputlist!(modelinputparameters,Parameter("K_Mg2PPi",10^(2.57),10^(1.57),10^(4.8),0.58))
    addtoinputlist!(modelinputparameters,Parameter("K_MgPi",10^(1.88)))
    addtoinputlist!(modelinputparameters,Parameter("Mg2PPi_eq",1.4e-5,1.4e-10,1.4e-2,0.3))

    return IVTmodel(modelinputparameters)
end