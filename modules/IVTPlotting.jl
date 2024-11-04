RNAindex = 2
PPiindex = 3
Mgindex = 8
alpha = 0.05

"""
    parseexceltuple(s::AbstractString) 

Parse tuple cells in csv data.
"""
function parseexceltuple(s::AbstractString) 
    if s[1]=='('
        if length(s)==2
            return []
        elseif s[end-1]==','
            return parse.(Float64, split(s[2:end-2], ','))
        else
            return (parse.(Float64, split(s[2:end-1], ',')))
        end
    else
          return [parse(Float64, s)]
    end
end
function parseexceltuple(s::Number) 
    return [s]
end

function parseexceltuple(s::Missing)
    return []
end

"""
   plotfromcsv(model,paramslist,parametercovariancematrix,filename; kwargs...)

Generate plot of model predictions alongside data. Uses csv file name.
# Keyword Arguments
- `labels::Array(String)`: labels to add to plot legend.
- `range::Array(Integer)`: rows of csv file to add to plot. Defaults to all rows.
- `plotsize = (400,400)`
- `multiplot = true`: Add a new plot window for each reaction condition.
- `maximumyield = true`: Add line denoting maximum yield of reaction.
- `multiplemaximum = false`: Give each maximum yield line a new color and legend label.
- `dataerrorbars = false`: Add error bars to data points from csv.
- `mcuncertainty = true`: Use Monte Carlo sampling from covariance matrix to generate prediction interval. If false uses linearized approximation.
- `nmc = 1000`: Number of Monte Carlo samples to use for uncertainty calculations.
- `plotsplit = (0,0)`: Start color scheme at point i of color gradient with maximum index j. Input as plotsplit = (i,j)
...
"""
function plotfromcsv(model,paramslist,parametercovariancematrix,filename; kwargs...)
    df = CSV.read(filename, DataFrame)
    rawdatamat = Matrix(df)
    return plotdata(model,paramslist,parametercovariancematrix,rawdatamat; kwargs...)
end


function plotdata(model,paramslist,parametercovariancematrix,data;labels = [], range = 1:size(data)[1], plotsize = (400,400), multiplot = true, maximumyield=true, multiplemaximum = false,plotsplit = (0,0),dataerrorbars = false,mcuncertainty = true, nmc = 1000, markersize = 3)
    params = fullparameterset(model,paramslist)
    
    #Multiplot creates a set of plot windows. Turning multiplot off plots everything in the same window.
    if multiplot
        pltvec = [plot() for i in 1:length(range)]
    else
        masterplot = plot()
    end
   
    init_speciesindex = data[1,1]
    if init_speciesindex == 2
        y_label = "RNA Yield (μM)"
    else
        y_label = "Concentration (mM)"
    end

    for (plotind,i) in enumerate(range)
        if multiplot
            plt = plot()
        end
        color = cgrad(:berlin, max(2,plotsplit[1]+length(range)), categorical = true)[plotind+plotsplit[2]]
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

        speciesindex = Int(data[i,1])
        times = parseexceltuple(data[i,2])
        yields = parseexceltuple(data[i,19])
        
        confinput = parseexceltuple(data[i,20])
        #Can either input uncertainty as scalar or vector
        if length(confinput) == 1
            confint = confinput .* ones(length(yields))
        else
            confint = confinput
        end

        #Adding error bars
        if dataerrorbars
            confint = confint
        else
            confint = 0
        end

        if speciesindex == 2
            speciesoutputfunction = sol -> sol[speciesindex,:]
            scalingfactor = 1e6
            y_label = "RNA Yield (μM)"
        elseif speciesindex != 11
            speciesoutputfunction = sol -> sol[speciesindex,:]
            scalingfactor = 1e3
            y_label = "Concentration (mM)"
        else
            speciesoutputfunction = sol -> -log10.(sol[speciesindex,:])
            scalingfactor = 1
            y_label = "pH"
        end
        

        stoich = SVector(N_A, N_U, N_C, N_G)

        inputs =(T7RNAP = T7RNAP, ATP = ATP,UTP = UTP,CTP = CTP,GTP = GTP, Mg = Mg, Buffer = Buffer, DNA = DNA, final_time = 1.2*maximum(times))
    
        sol = runDAE_batch(params, inputs, PPiase = PPiase, stoich = stoich, Cap = Cap)
        plot!(sol.t,(speciesoutputfunction(sol)) .* scalingfactor,color = color,label = "",linewidth =3)

        tvals = sol.t
        maximumRNA = scalingfactor*min((ATP/N_A),(UTP/N_U),(CTP/N_C),(GTP/N_G))
        predictionerror = zeros(length(tvals))


        #Generate uncertainty in predictions
        if mcuncertainty # Using mc sampling for uncertainty
            mcensemble = zeros(nmc,length(tvals))
            mean = paramslist
            d = MvNormal(mean, Hermitian(parametercovariancematrix))
            for i in 1:nmc
                x = rand(d, 1)
                sampleparams = fullparameterset(model,x)
                sol = runDAE_batch(sampleparams, inputs, PPiase = PPiase, stoich = stoich, Cap = Cap)
                timepointsol = sol(tvals)
                mcensemble[i,:] = (speciesoutputfunction(timepointsol)) .* scalingfactor
            end
            lower_pointwise_CB = [percentile(mcensemble[:,j],100*alpha/2) for j in 1:length(tvals)]
            upper_pointwise_CB = [percentile(mcensemble[:,j],100*(1-alpha/2)) for j in 1:length(tvals)]

        else #Using Linearized approximation for uncertainty calculation
            predictionerror = predictionuncertainty_linearized_multipoint(speciesindex,alpha,model,params,parametercovariancematrix, inputs, tvals)
            lower_pointwise_CB = max.(scalingfactor*(sol[speciesindex,:]) .- scalingfactor*predictionerror,zeros(length(sol.t)))[:,:]
            upper_pointwise_CB = min.(scalingfactor*(sol[speciesindex,:]) .+ scalingfactor*predictionerror,maximumRNA)[:,:]
        end
        
        plot!(tvals, lower_pointwise_CB, fillrange = upper_pointwise_CB, fillalpha = 0.35,alpha=0.0, color = color,linewidth = 0.0,label="")
        
        #Shows purple line to denote maximum yield
        if maximumyield 
            if multiplot || !multiplemaximum
                plotlabel = ""
                plotcolor = palette(:tab10)[5]
            else
                plotlabel = "Maximum Yield: "
                plotcolor = color
            end
            plot!(sol.t,maximumRNA .* ones(length(sol.t)),color = plotcolor,linewidth = 4,linestyle = :dash,label = plotlabel)
        end
        
        #Uses labels input for plot legend
        if isempty(labels)
            label = "Mg:"*string(Int(round(1000*Mg)))*" mM"
        else
            label = labels[plotind]
        end

        scatter!(times,yields,mc = color,markersize = markersize,label = label, yerror = confint,markerstrokecolor=palette(:grays,10)[3])
        if multiplot
            pltvec[plotind] = plt
        end
    end
    if multiplot
        if length(pltvec)<4
            masterplot = plot(pltvec...,size = plotsize,xtickfontsize=12,ytickfontsize=12,xguidefontsize=15,yguidefontsize=15,grid = false,layout = (1,3),legend=:outerright,xlabel = "Time (h)",ylabel = y_label)
        else
            masterplot = plot(pltvec...,size = plotsize,xtickfontsize=12,ytickfontsize=12,xguidefontsize=15,yguidefontsize=15,grid = false,legend=:outerright,xlabel = "Time (h)",ylabel = y_label)
        end
    else
        plot!(size = plotsize,xtickfontsize=12,ytickfontsize=12,xguidefontsize=15,yguidefontsize=15,grid = false,legend=:outerright,xlabel = "Time (h)",ylabel = y_label)
    end
    plot!(ylims=(0,Inf))
    return masterplot
end
