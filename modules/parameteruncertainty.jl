"""
    filtersensitivitymatrix(u,mask)

Take sensitivity matrix and mask, filter matrix to only contain rows representing fitted parameters.
""" 
function filtersensitivitymatrix(sensitivitymat,mask)
    result = zeros(Int(sum(mask)),size(sensitivitymat)[2])
    sumcounter = 1
    for i in 1:size(sensitivitymat)[1]
        if mask[i]==1
            result[sumcounter,:] = sensitivitymat[i,:]
            sumcounter+=1
        end
    end
    result
end

"""
    getcorrelationmatrix(mat)

Take covariance matrix and return correlation matrix.
""" 
function getcorrelationmatrix(mat)
    D = sqrt.(diagm(diag(mat)))
    DInv = inv(D);
    return DInv * mat * DInv
end

"""
    parametermask(model)

Take model, return prior component of information matrix for use in covariance matrix calculation.
""" 
function priormatrix(model)
    mask = parametermask(model)
    priorlen = Int(sum(mask))
    priormatrix = zeros(priorlen,priorlen)
    sumval = 1
    for i in 1:length(model.parameters)
        if model.parameters[i].hasprior
            priormatrix[sumval,sumval] = 1/model.parameters[i].variance
        end
        if model.parameters[i].isfitted
            sumval+=1
        end

    end
    return priormatrix
end

"""
    getcovariancematrix(model,data,optparamlist)

Take list of parameters from optimizer and return covariance matrix. 
""" 
function getcovariancematrix(model,data,optparamlist; customfile = false,customfilename = "")
    mask = parametermask(model)
    u,i = sensitivitymatricies(model,data,optparamlist,false)#Concentation trajectories of RNA and PPi
    ini = initialratesensitivity(model,data,optparamlist,false)#Initial transcription rate data
    tit = PPititrationsensitivity(model,data,optparamlist,false)#Mg2PPi solubility data

    filtered_u = filtersensitivitymatrix(u,mask)
    filtered_i = filtersensitivitymatrix(i,mask)
    filtered_ini = filtersensitivitymatrix(ini,mask)
    filtered_tit = filtersensitivitymatrix(tit,mask)

    u_covcomp = filtered_u*diagm(1 ./(reshape(data.RNAstdev',(13*9))) .^2)*filtered_u'
    i_covcomp = filtered_i*diagm(1 ./(reshape(data.PPistdev',(13*9))) .^2)*filtered_i'
    titration_covcomp = filtered_tit*inv(diagm(data.Mgstdev.^2))*filtered_tit'
    initial_covcomp = filtered_ini*inv(diagm(data.initialrateyieldstdev.^2))*filtered_ini'
    prior_covcomp = priormatrix(model).^2

    if customfile
        df = CSV.read(customfilename, DataFrame)
        data = Matrix(df)
        custom_sm, stdevvector = customsensitivitymatrix(model,data,optparamlist,false)
        custom_sm_filtered = filtersensitivitymatrix(custom_sm,mask)
        custom_covcomp = custom_sm_filtered*inv(diagm(stdevvector .^2))*custom_sm_filtered'
    else
        custom_covcomp = zeros(Int(sum(mask)),Int(sum(mask)))
    end
    parametercovariancematrix = inv(u_covcomp .+ i_covcomp.+ initial_covcomp .+ prior_covcomp .+ titration_covcomp .+ custom_covcomp)
end

"""
    samplecovariancematrix(covariancematrix)

Calculate 2sigma confidence intervals of each parameter by sampling covariance matrix.
""" 
function samplecovariancematrix(covariancematrix)
    mean = zeros(size(covariancematrix)[1])
    d = MvNormal(mean, Hermitian(covariancematrix))
    x = rand(d, 10000000)

    samplingresults = zeros(size(covariancematrix)[1])
    for i in 1:length(samplingresults)
        samplingresults[i] = Statistics.var(x[i,:])
    end
    return sqrt.(samplingresults)*2
end