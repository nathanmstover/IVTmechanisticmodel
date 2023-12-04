"""
    localopt(model,data,time,start,residualevaluation;verbose = false)

Perform local optimization of model from start point based on residualevaluation function.
"""
function localopt(model,data,time,start,residualevaluation;verbose = false)
    origin, lower_bounds, upper_bounds = generateoptimizationsettings(model)
    num = length(start)

    opt = Opt(:LD_LBFGS, num)   
    opt.ftol_rel = 1e-6
    opt.min_objective = (x,g)->gradoptimizationfunctionwrapper(x,g,model,data,residualevaluation;verbose = verbose)

    opt.lower_bounds = lower_bounds
    opt.upper_bounds = upper_bounds

    opt.maxtime = time# seconds

    (minf,minx_log,ret) = NLopt.optimize(opt, start)
    minx = [10^i for i in minx_log]
    numevals = opt.numevals # the number of function evaluations
    if verbose
        println("got $minf at $minx after $numevals iterations (returned $ret)")
    end
    return numevals,minf,start,minx_log,Int64(String(ret)[1])
end

"""
    gradoptimizationfunctionwrapper(fititerationparameters::Vector, grad::Vector, model, data,residualevaluation;verbose = false)

Evaluate objective function and gradients for optimization algorithm.
"""
function gradoptimizationfunctionwrapper(fititerationparameters::Vector, grad::Vector, model, data,residualevaluation;verbose = false)

    func = x -> residualevaluation(model,data,x)
    result = DiffResults.GradientResult(fititerationparameters)
    result = ForwardDiff.gradient!(result, func, fititerationparameters)
    autodiffgrad = DiffResults.gradient(result)
    for i in 1:(length(grad))
        grad[i] = autodiffgrad[i]
    end

    residual = DiffResults.value(result)
    if true#verbose
        println(residual)
    end
    return residual 
end

"""
    generateoptimizationsettings(model)

Generate start point for optimization based on initial guesses entered in model construction."""
function generateoptimizationsettings(model)
    start::Vector{Float64} = []
    lowerbound::Vector{Float64} = []
    upperbound::Vector{Float64} = []
    for i in model.parameters
        if i.isfitted
                append!(start,[log10(i.value)])  
                append!(lowerbound,[log10(i.upperbound)])  
                append!(upperbound,[log10(i.lowerbound)])  
        end
    end
    return start,lowerbound,upperbound
end
