using Pkg
Pkg.activate("IVTmodel")
Pkg.instantiate()
using DifferentialEquations
using Plots
using NLopt
using LinearAlgebra
using Metrics
using CSV, DataFrames
using StaticArrays
using ForwardDiff
using NamedTupleTools
using BenchmarkTools
using ComponentArrays
using Distributions
using Statistics
using GenericLinearAlgebra
using DiffResults
using FLoops
using SharedArrays
using Plots.PlotMeasures
using FlexiMaps
using LaTeXStrings
using DelimitedFiles
using StatsBase
using NLsolve

include("./model.jl")
include("./parameterbackend.jl")
include("./akamadataprocessing.jl")
include("./akamaplotting.jl")
include("./IVTplotting.jl")
include("./modelfitting.jl")
include("./residualcalculation.jl")
include("./parameteruncertainty.jl")
include("./predictionuncertainty.jl")
include("./akamasensitivity.jl")
include("./solubility.jl")


#Generates settings for parameters used
fittingmodel = setupmodel()
akamafittedparametersmatrix = Matrix(CSV.read("outputs/fittedparameters.csv", DataFrame,header=false))
akamafittedparameterslist = reshape(akamafittedparametersmatrix,(size(akamafittedparametersmatrix)[1],))
covariancematrix = Matrix(CSV.read("outputs/covariancematrix.csv", DataFrame,header=false))
akamafittedparams = fullparameterset(fittingmodel,akamafittedparameterslist);
alpha = 0.05

"""
    parameterized_IVT_PFR(parameters,T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)

Takes reactor inputs, (T7RNAP, etc) along with parameter values to generate single prediction. Can represent batch reaction (where spacetime=reaction time) or idealized PFR (where spacetime= reactor length/fluid speed)
""" 
function parameterized_IVT_PFR(parameters,T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)
    inputs = (T7RNAP = T7RNAP, ATP = ATP,UTP = UTP,CTP = CTP,GTP = GTP, Mg = Mg, Buffer = 0.040, DNA = DNA, final_time = spacetime)
    sol = runDAE_batch(parameters, inputs; kwargs...)
    RNAoutput = 1e6*sol[2,end]
    initialNTP = (ATP+UTP+CTP+GTP)
    finalNTP = sum(sol[4:7,end])
    NTPconversion = 1 - finalNTP /initialNTP
    return [RNAoutput, NTPconversion]
end

"""
    IVT_PFR(T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)

Part of IVT API. Takes reactor inputs to generate single prediction using maximum likelihood parameter set. Can represent batch reaction (where spacetime=reaction time) or idealized PFR (where spacetime = reactor length/fluid speed)

### Inputs:
- T7RNAP: T7RNA polymerase input concentration (mol/L)
- ATP-GTP: Concentration of respective NTP in mol/L
- Mg: Concentration of Magnesium salt in mol/L
- DNA: Concentration of DNA in nanomoles/L
- spacetime: Time in hours

### Important keyword arguments:
- stoich = SVector(231, 246, 189, 202): tuple of integers describing number of A,U,C,G in target sequence.
- PPiase = 0.0: concentration of pyrophosphatase enzyme in U/uL
- Cap = 0.0: concentration of cap analogue in mol/L
- tol = 1e-5: float of DAE solver tolerance. 

### Outputs: 
- Vector of two values [RNAoutput, NTPconversion]
- RNAoutput = Concentration of RNA product in outlet stream in umol/L
- NTPconversion = fractional conversion of NTP feedstock (bounded between 0-1)

""" 
function IVT_PFR(T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)
    parameterized_IVT_PFR(akamafittedparams,T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)
end


"""
    IVT_PFR_uncertainty(T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; n_mc_samples = 1000, kwargs...)

Part of IVT API. Takes reactor inputs to generate median and 95% confidence bounds using maximum likelihood parameter set. Can represent batch reaction (where spacetime=reaction time) or idealized PFR (where spacetime = reactor length/fluid speed)

### Inputs:
- T7RNAP: T7RNA polymerase input concentration (mol/L)
- ATP-GTP: Concentration of respective NTP in mol/L
- Mg: Concentration of Magnesium salt in mol/L
- DNA: Concentration of DNA in nanomoles/L
- spacetime: Time in hours

### Important keyword arguments:
- stoich = SVector(231, 246, 189, 202): tuple of integers describing number of A,U,C,G in target sequence.
- PPiase = 0.0: concentration of pyrophosphatase enzyme in U/uL
- Cap = 0.0: concentration of cap analogue in mol/L
- tol = 1e-5: float of DAE solver tolerance. 

### Outputs: 
- tuple of vectors (maximumlikelihood, lowerCB, upperCB), representing uncertainty of output (median, lower 95% confidence bound, upper 95% confidence bound.) 
- Each vector comprises two values [RNAoutput, NTPconversion]
- RNAoutput = Concentration of RNA product in outlet stream in umol/L
- NTPconversion = fractional conversion of NTP feedstock (bounded between 0-1)


""" 
function IVT_PFR_uncertainty(T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; n_mc_samples = 1000, kwargs...)
    n_outputs = 2
    maximumliklihood = IVT_PFR(T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)
    mcensemble = zeros(2,n_mc_samples)
    mean = akamafittedparameterslist
    d = MvNormal(mean, Hermitian(covariancematrix))
    for i in 1:n_mc_samples
        x = rand(d, 1)
        sampleparams = fullparameterset(fittingmodel,x)
        sol = parameterized_IVT_PFR(sampleparams,T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)
        mcensemble[:,i] = sol
    end
    maximumliklihood = IVT_PFR(T7RNAP,ATP,UTP,CTP,GTP,Mg,DNA,spacetime; kwargs...)
    lowerCB = [percentile(mcensemble[i,:],100*alpha/2) for i in 1:n_outputs]
    upperCB = [percentile(mcensemble[i,:],100*(1-alpha/2)) for i in 1:n_outputs]
    return maximumliklihood, lowerCB, upperCB
end

