#Global variables for species indexing
RNAindex = 2
PPiindex = 3
Mgindex = 8
#Global variable defines 95% prediction intervals (99% -> alpha = 0.01)
alpha = 0.05

"""
    plotinitialrates(model,data,params, parametercovariancematrix, showconfidence)  

Generate figure 3A of main text showing model fit to initial rate Akama data."""
function plotinitialrates(model,data,params, parametercovariancematrix, showconfidence)
    plt = plot(legend_title = "Mg",legend_title_font_pointsize = 9)
    NTPpoints = 100
    NTPinput = LinRange(1e-18,0.02,NTPpoints)
    mask = parametermask(model)
    for (Mgindex,Mgconc) in enumerate([0.002,0.004,0.008,0.014,0.020])
        color = cgrad(:berlin, 5, categorical = true)[Mgindex]
        
        RNAyield = zeros(NTPpoints)
        predictionerror = zeros(NTPpoints)
        for (NTPindex,NTPconc) in enumerate(NTPinput)
            inputs =(T7RNAP = 1e-7, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = Mgconc, Buffer = 0.040, DNA = 7.4, final_time = 0.08333333333)

            sol = runDAE_batch(params,inputs)
            RNAyield[NTPindex] = sol.u[end][RNAindex]
            if showconfidence
                predictionerror[NTPindex] = predictionuncertainty_linearized(RNAindex,alpha,model,params,parametercovariancematrix, inputs; saveevery = false)
            else
                predictionerror[NTPindex] = 0
            end
        end
        plot!(1000*NTPinput,1e6*RNAyield,linewidth = 2.5,color = color, label = "",z_order = :back)

        lower_pointwise_CB = max.(1e6*RNAyield .- 1e6*predictionerror,0)
        upper_pointwise_CB = 1e6*RNAyield .+ 1e6*predictionerror
        plot!(1000*NTPinput, lower_pointwise_CB[:,:], fillrange = upper_pointwise_CB[:,:], fillalpha = 0.35,alpha=0.0, color = color,linewidth = 0.0,label="",z_order = :back)
       
        
        exptdatainput = data.initialrateconcentrationinputs[(data.initialrateconcentrationinputs[:,2] .== Mgconc),:]
        exptdataoutput = data.initialrateyieldoutputs[(data.initialrateconcentrationinputs[:,2] .== Mgconc)]
        exptdatastdev = data.initialrateyieldstdev[(data.initialrateconcentrationinputs[:,2] .== Mgconc)]

        scatter!(1000*exptdatainput[:,1],1e6*exptdataoutput,yerr = 1e6*exptdatastdev,markersize = 4,palette =:Dark2_5,label=string(Int(round(Mgconc*1000, sigdigits = 2)))*" mM",mc = color)
    end
    plot!(xlabel = "NTP (mM)",ylabel = "RNA (μM)",linewidth=3.5,legend=:outerright,ylims = (0,Inf))
    return plt   
end

"""
    function plotMgtitrations(model,data,params, parametercovariancematrix, showconfidence, Mg = 0.004, time = 24, showdata = true, maxNTP = 0.8)

Generate figure 3B of main text showing model fit to Mg2PPi solubility Akama data."""
function plotMgtitrations(model,data,params, parametercovariancematrix, showconfidence, Mg = 0.004, time = 24, showdata = true, maxNTP = 0.8)
    PPipoints = 100
    addedNTP = (Mg/0.004) * 4 .* LinRange(1e-7,maxNTP*1e-3,5)
    addedPPi = LinRange(1e-18,(Mg/0.004) * 0.009,PPipoints)
    mask = parametermask(model)
    plt = plot(xlabel = "PPi (mM)", ylabel = "Mg (mM)", legend_title = "NTP",legend_title_font_pointsize = 9)
    for (NTPindex,NTPconc) in enumerate(addedNTP)
        color = cgrad(:berlin, 5, categorical = true)[NTPindex]
        totalMgoutput = zeros(PPipoints)
        predictionerror = zeros(PPipoints)
        for (PPiind,PPiconc) in enumerate(addedPPi)
            inputs =(T7RNAP = 0, ATP = NTPconc/4,UTP = NTPconc/4,CTP = NTPconc/4,GTP = NTPconc/4, Mg = Mg, Buffer = 0.040, DNA = 1e-8, final_time = time)
            sol = runDAE_batch(params,inputs; PPi = PPiconc, saveevery = false)
            totalMgoutput[PPiind] = sol.u[end][Mgindex]
            if showconfidence
                predictionerror[PPiind] = predictionuncertainty_linearized(Mgindex,alpha,model,params,parametercovariancematrix, inputs; PPi = PPiconc, saveevery = false)
            else
                predictionerror[PPiind] = 0
            end
        end
        plot!(1000*addedPPi,1000*totalMgoutput,linewidth = 2.5,color = color,label = string(round(NTPconc*1000,digits = 2))*" mM",z_order = :back)
        
        lower_pointwise_CB = max.(1000*totalMgoutput .- 1000*predictionerror,0)
        upper_pointwise_CB = 1000*totalMgoutput .+ 1000*predictionerror
        plot!(1000*addedPPi, lower_pointwise_CB[:,:], fillrange = upper_pointwise_CB[:,:], fillalpha = 0.35,alpha=0.0, color = color,linewidth = 0.0,label="",z_order = :back)
        
        if showdata
            PPidata = data.Mgconcentrationinputs[1:10,2]
            Mgdata = data.Mgoutputs[(10*(NTPindex-1)+1):(10*(NTPindex-1)+10),1]
            Mgstdevdata = data.Mgstdev[(10*(NTPindex-1)+1):(10*(NTPindex-1)+10),1]
            scatter!(1000*PPidata, 1000*Mgdata, yerror = 1000*Mgstdevdata, markersize = 4,mc = color,label = "")
        end
    end
    return plt
end

"""
    function gendatasubplot(parametervaried,output,model, data, parameters, parametercovariancematrix, genmodelpredictions,showconfidence)

Generate individual windows of figure 2 in the main text showing model fit to Mg2PPi solubility Akama data.""" 
function gendatasubplot(parametervaried,output,model, data, parameters, parametercovariancematrix, genmodelpredictions,showconfidence)
    mask = parametermask(model)
    parameters = ComponentArray(parameters)
    plt = plot()
    maximumRNA = 1e6*(0.0008/246)#1e6*(0.0032/868)
    
    if parametervaried == "RNAP"
        paramcolumn = 1
    elseif parametervaried == "NTP"
        paramcolumn = 2
    elseif parametervaried == "Mg"
        paramcolumn = 3
    else
        throw(DomainError(parametervaried, "Parameter to vary not listed"))
    end
    
    if output == "RNA"
        dataind = 5
        stdind = 6
        modelind = RNAindex
        y_label = "RNA (μM)"
        scale_factor = 1e6
        ylims = (0,maximumRNA+0.1)
    elseif output == "PPi"
        dataind = 7
        stdind = 8
        modelind = PPiindex
        scale_factor = 1e3
        y_label = "PPi (mM)"
        ylims = (0,3)

    else 
        throw(DomainError(parametervaried, "Parameter to vary not listed"))
    end
    
    
    for (colorind,k) in enumerate(sort(unique(data.concentrationinputs[:,paramcolumn])))
        color = cgrad(:berlin, 5, categorical = true)[colorind]
        
        if (parametervaried == "RNAP" && k == 0.1*1e-6)||(parametervaried == "NTP" && k == 0.0032)||(parametervaried == "Mg" && k == 0.008)
            indexbitvec = [(data.concentrationinputs[:,1] .== 0.1*1e-6) .& (data.concentrationinputs[:,2].==0.0032) .& (data.concentrationinputs[:,3].==0.008),:][1]
            indexval = dot(indexbitvec,collect(1:13))
        else
            indexbitvec = [(data.concentrationinputs[:,paramcolumn] .== k),:][1]
            indexval = dot(indexbitvec,collect(1:13))
        end
        if output == "PPi"
            plot!(legend_title_font_pointsize = 9)
            if parametervaried == "RNAP"
                legend_label = [string(round(i*1e6,sigdigits = 2))*" μM" for i in k]
                legend_title = "T7 RNAP"
            elseif parametervaried == "Mg"
                legend_label = [string(Int(round(i*1000)))*" mM" for i in k]    
                legend_title = "Mg"

            else
                legend_label = [string(round(i*1000,sigdigits = 2))*" mM" for i in k]
                legend_title = "NTP (total)"
            end
        else
            legend_label = ["" for i in k]
        end

        if output == "RNA"
            if genmodelpredictions
                scatter!(data.timeinputs[indexval,:],1e6*data.RNAyieldoutputs[indexval,:],yerror=1e6*data.RNAstdev[indexval,:],markersize = 4,mc=color,label=legend_label,legend_column = 1,markerstrokewidth = 0.3)#, title = "Varying "*parametervaried*": "*output)
            else
                plot!(data.timeinputs[indexval,:],1e6*data.RNAyieldoutputs[indexval,:],yerror=1e6*data.RNAstdev[indexval,:],markersize = 4,mc=color,label=legend_label,lc=color,linewidth=3.5,z_order = :back)#, title = "Varying "*parametervaried*": "*output)
            end
        else
            if genmodelpredictions
                scatter!(data.timeinputs[indexval,:],1e3*data.PPiyieldoutputs[indexval,:],yerror=1e3*data.PPistdev[indexval,:],markersize = 4,mc=color,label=legend_label,legend_title=legend_title, legend_column = 1,markerstrokewidth = 0.3)#, title = "Varying "*parametervaried*": "*output)
            else
                plot!(data.timeinputs[indexval,:],1e3*data.PPiyieldoutputs[indexval,:],yerror=1e3*data.PPistdev[indexval,:],markersize = 4,mc=color,label=legend_label,lc=color,linewidth=3.5,z_order = :back)#, title = "Varying "*parametervaried*": "*output)
            end
        end
            
        if genmodelpredictions
            RNAP = 0.1*1e-6
            NTP = 0.0032
            Mg = 0.008
            if parametervaried == "RNAP"
                RNAP = k
            elseif parametervaried == "NTP"
                NTP = k
            elseif parametervaried == "Mg"
                Mg = k
            else
                throw(DomainError(parametervaried, "Parameter to vary not listed"))
            end

            inputs =(T7RNAP = RNAP, ATP = NTP/4,UTP = NTP/4,CTP = NTP/4,GTP = NTP/4, Mg = Mg, Buffer = 0.040, DNA = 7.4, final_time = 2)
            sol = runDAE_batch(parameters,inputs)

            plot!(sol.t,scale_factor*sol[modelind,:],lc=color,label="",linewidth=2.5,z_order = :back)
            
            tvals = sol.t

            if showconfidence
                if output == "RNA"
                    predictionerror = predictionuncertainty_linearized_multipoint(RNAindex,alpha,model,parameters,parametercovariancematrix, inputs, tvals)
                else
                    predictionerror = predictionuncertainty_linearized_multipoint(PPiindex,alpha,model,parameters,parametercovariancematrix, inputs, tvals)
                end
            else
                predictionerror = zeros(length(tvals))
            end

            lower_pointwise_CB = max.(scale_factor*sol[modelind,:] .- scale_factor*predictionerror,zeros(length(sol.t)))
            upper_pointwise_CB = scale_factor*sol[modelind,:] .+ scale_factor*predictionerror

            plot!(sol.t, lower_pointwise_CB[:,:], fillrange = upper_pointwise_CB[:,:], fillalpha = 0.35,alpha=0.0, color = color,linewidth = 0.0,label="",z_order = :back)
        end
    end
    plot!(legend=:outertop,ylim = ylims)
    return plt

end

"""
    function revision1plot1(model, data,parameters, parametercovariancematrix, genpredictions,showconfidence)

Generate figure 2 of main text.""" 
function revision1plot1(model, data,parameters, parametercovariancematrix, genpredictions,showconfidence)
    plt = plot()
    widthfactor = 1
    l = @layout([grid(1, 3){0.575h}; grid(1, 3)])
    maximumRNA = 1e6*(0.0008/246)#1e6*(0.0032/868)
    p1 = gendatasubplot("RNAP","PPi",model,data, parameters, parametercovariancematrix, genpredictions,showconfidence)
    plot!(ylabel = "PPi (mM)",leftmargin = 5mm)
    p2 = gendatasubplot("NTP","PPi",model,data, parameters, parametercovariancematrix, genpredictions,showconfidence)
    p3 = gendatasubplot("Mg","PPi",model,data, parameters, parametercovariancematrix, genpredictions,showconfidence)

    p4 = gendatasubplot("RNAP","RNA",model, data,parameters, parametercovariancematrix, genpredictions,showconfidence)
    plot!(LinRange(0,2,100),ones(100) .* maximumRNA,color = palette(:tab10)[5],linewidth = 4,linestyle = :dash,label = "")
    plot!(ylabel = "RNA (μM)",leftmargin = 5mm,bottommargin = 10mm)
    plot!(xlabel = "Time (hours)")
    p5 = gendatasubplot("NTP","RNA",model,data, parameters, parametercovariancematrix,genpredictions,showconfidence)
    plot!(bottommargin = 10mm)
    plot!(xlabel = "Time (hours)")
    p6 = gendatasubplot("Mg","RNA",model,data, parameters, parametercovariancematrix, genpredictions,showconfidence)
    plot!(bottommargin = 10mm)
    plot!(xlabel = "Time (hours)")
    plot!(LinRange(0,2,100),ones(100) .* maximumRNA,color = palette(:tab10)[5],linewidth = 4,linestyle = :dash,label = "")


    ex_title_left = Plots.plot(title = "Left group title", grid = false, showaxis = false)
    ex_title_right = Plots.plot(title = "Right group title", grid = false, showaxis = false)
    
    plt1 = plot!(p1, p2, p3, p4, p5, p6, layout = l, size = (950,800))
    plot!(xtickfontsize=12,ytickfontsize=12,xguidefontsize=13,yguidefontsize=13,grid = false)
    fontsize=24
    
    annotate!(sp=1,[(relativex(-0.07; sp=1), relativey(1.11; sp=1), text("A",fontsize))])
    annotate!(sp=2,[(relativex(-0.07; sp=2), relativey(1.11; sp=2), text("B",fontsize))])
                                                        
    annotate!(sp=3,[(relativex(-0.07; sp=3), relativey(1.11; sp=3), text("C",fontsize))])
    annotate!(sp=4,[(relativex(-0.07; sp=4), relativey(1.08; sp=4), text("D",fontsize))])
                                                    
    annotate!(sp=5,[(relativex(-0.07; sp=5), relativey(1.08; sp=5), text("E",fontsize))])
    annotate!(sp=6,[(relativex(-0.07; sp=6), relativey(1.08; sp=6), text("F",fontsize))])
                                                    
    return plt1
end


"""
    function revision1plot2(model, data,parameters, parametercovariancematrix, genpredictions,showconfidence)

Generate figure 3 of main text.""" 
function revision1plot2(model, data,parameters, parametercovariancematrix, genpredictions,showconfidence)
    plt = plot()
    widthfactor = 1
    #l = @layout([[A B C; D E F] [G; H]])grid(4, 2){0.33w}
    l = @layout([A B])
    maximumRNA = 1e6*(0.0008/246)#1e6*(0.0032/868)

    p2 = plotMgtitrations(model,data, parameters,parametercovariancematrix, showconfidence)
    plot!(bottommargin = 5mm,topmargin = 5mm)
    p1 = plotinitialrates(model,data,parameters,parametercovariancematrix, showconfidence)
    plot!(bottommargin = 9mm)

    ex_title_left = Plots.plot(title = "Left group title", grid = false, showaxis = false)
    ex_title_right = Plots.plot(title = "Right group title", grid = false, showaxis = false)
    
    plt1 = plot!(p1, p2, layout = l, size = (800,400),leftmargin = 5mm,topmargin = 8mm)
    plot!(xtickfontsize=12,ytickfontsize=12,xguidefontsize=13,yguidefontsize=13,grid = false)
    fontsize=24
    


    
    annotate!(sp=1,[(relativex(-0.13; sp=1), relativey(1.08; sp=1), text("A",fontsize))])
    annotate!(sp=2,[(relativex(-0.13; sp=2), relativey(1.08; sp=2), text("B",fontsize))])

    return plt1
end

"""
    relative(f, r; sp)

Generate relative coordinates for plotting.""" 
function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2]-lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)

"""
    plotmodelpredictionPPiase(params, parametercovariancematrix,ppimodel,PPiase,model)

Generate model predictions for main text figure 4.""" 
function plotmodelpredictionPPiase(params, parametercovariancematrix,ppimodel,PPiase,model)
    if ppimodel
        colorincRNA = 4
        coloringPPi = 3
    else
        colorincRNA = 1
        coloringPPi = 2
    end
    color = palette(:RdYlGn_4)

    inputs =(T7RNAP = 1e-7, ATP = 0.0008,UTP = 0.0008,CTP = 0.0008,GTP = 0.0008, Mg = 0.008, Buffer = 0.040, DNA = 7.4, final_time = 1.0)
    sol = runDAE_batch(params, inputs, PPiase = PPiase)
    tvals = sol.t

    RNAmodelpredictionerror = predictionuncertainty_linearized_multipoint(RNAindex,alpha,model,params,parametercovariancematrix, inputs, tvals; PPiase)
    PPimodelpredictionerror = predictionuncertainty_linearized_multipoint(PPiindex,alpha,model,params,parametercovariancematrix, inputs, tvals; PPiase)

    RNAlower_pointwise_CB = 1e6*sol[RNAindex,:] .- 1e6*RNAmodelpredictionerror
    RNAupper_pointwise_CB = 1e6*sol[RNAindex,:] .+ 1e6*RNAmodelpredictionerror

    PPilower_pointwise_CB = 1e3*sol[PPiindex,:] .- 1e3*PPimodelpredictionerror
    PPiupper_pointwise_CB = 1e3*sol[PPiindex,:] .+ 1e3*PPimodelpredictionerror


    plot!(sol.t,1e6*sol[RNAindex,:],label="",linewidth=3,color = color[colorincRNA])
    plot!(sol.t, RNAlower_pointwise_CB[:,:], fillrange = RNAupper_pointwise_CB[:,:], fillalpha = 0.35,alpha=0.0, color = color[colorincRNA],linewidth = 0.0,label="")

    plot!(sol.t,1e3*sol[PPiindex,:],label="",linewidth=3,color = color[coloringPPi])
    plot!(sol.t, PPilower_pointwise_CB[:,:], fillrange = PPiupper_pointwise_CB[:,:], fillalpha = 0.35,alpha=0.0, color = color[coloringPPi],linewidth = 0.0,label="")
end

"""
    plotakamaPPiase(PPiasedata,params, parametercovariancematrix,model)

Generate main text figure 4.""" 
function plotakamaPPiase(PPiasedata,params, parametercovariancematrix,model)
    PPiplotcolors = palette(:RdYlGn_4)
    plt = plot()
    maximumRNA = 1e6*(0.0008/246)
    plot!(LinRange(0,1,100),ones(100) .* maximumRNA,color = palette(:tab10)[5],linewidth = 4.5,linestyle = :dash,label = "Maximum Yield")
    scalingfactor = 1.6
    scatter!(PPiasedata[:,1],scalingfactor*1e6*PPiasedata[:,2],yerror=(4.3/sqrt(3))*1.6*1e6*PPiasedata[:,3],markersize = 4,label="RNA without PPiase",color = PPiplotcolors[1])
    scatter!(PPiasedata[:,1],1e3*PPiasedata[:,4],yerror=(4.3/sqrt(3))*1e3*PPiasedata[:,5],markersize = 4,label="PPi without PPiase",color = PPiplotcolors[2])
    plotmodelpredictionPPiase(params, parametercovariancematrix,false,0.0,model)
    scatter!(PPiasedata[:,1],scalingfactor*1e6*PPiasedata[:,6],yerror=(4.3/sqrt(3))*1.6*1e6*PPiasedata[:,7],markersize = 4,label="RNA with PPiase",color = PPiplotcolors[4])
    scatter!(PPiasedata[:,1],1e3*PPiasedata[:,8],yerror=(4.3/sqrt(3))*1e3*PPiasedata[:,9],markersize = 4,label="PPi with PPiase",color = PPiplotcolors[3])
    plotmodelpredictionPPiase(params, parametercovariancematrix,true,0.07,model)


    plt = plot!(xtickfontsize=12,ytickfontsize=12,xguidefontsize=15,yguidefontsize=15,grid = false,legend = :outerright)

    plot!(xlabel = "Time (hours)",ylabel = "RNA (μM), PPi, (mM)",size = (550,350),ylim = (0,maximumRNA+0.5))
    savefig(plt,"figures/highlight2.png")
    plt
end
