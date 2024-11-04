
"""
    initialconcentrationresidual!(F, x, ptot, param)

Generate residual for nonlinear solving of free ion concentrations. Only used for initialization of DAE.

"""
function initialconcentrationresidual!(F, x, ptot, param)
    H = exp(x[1])
    Mg = exp(x[2])

    (ions,(Chargebalance,Mgbalance)) = speciationmodel(param, 0, ptot.NTPtot, ptot.Mgtot, ptot.Buffertot, ptot.PPitot, ptot.Pitot, H, Mg, ptot.Na, ptot.Cl, ptot.OAc)
    #Algebraic Equations
    # H mass balance
    F[1] = Chargebalance
    # Mg mass balance
    F[2] = Mgbalance
    nothing
end

"""
    getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot)  where {T<:Real}

Take parameters and total concentration of ions, return free concentrations of H and Mg. Performs solving of nonlinear system of equations in log space. Only used for initialization of DAE.

"""
function getfreeconcentrations(params::AbstractArray{T}, NTPtot, Mgtot, Buffertot, PPitot, Pitot, Na, Cl, OAc)  where {T<:Real}
    #Defining Parameters
    initialtotalconcentrations = (NTPtot=NTPtot, Mgtot=Mgtot, Buffertot=Buffertot, PPitot = PPitot, Pitot = Pitot, Na = Na, Cl = Cl, OAc = OAc)
    guessfreeconcentrations = zeros(T,2)
    guessfreeconcentrations += log.([1.32e-8,0.00397])
    logsolvedfreeconcentrations = nlsolve((F,x)->initialconcentrationresidual!(F, x, initialtotalconcentrations, params), guessfreeconcentrations,ftol = 1e-8)
    solvedfreeconcentrations = [exp(x) for x in logsolvedfreeconcentrations.zero]
    return solvedfreeconcentrations
end

"""
    speciationmodel(param, NTPtot, Mgtot, Buffertot, PPitot, Pitot, H, Mg)

Take parameters and total concentration of ions, return concentrations of ionic species and ion equilibrium residuals for use in rate calculations and DAE solving.

"""
function speciationmodel(param, RNAtotalbases, NTPtot, Mgtot, Buffertot, PPitot, Pitot, H, Mg, Na, Cl, OAc;buffer_pka = 8.1)

    #Dimensionless concentrations for NTP
    HNTP_dimless = H*param.K_HNTP
    HMgNTP_dimless = HNTP_dimless*Mg*param.K_HMgNTP
    MgNTP_dimless = Mg*param.K_MgNTP
    Mg2NTP_dimless = MgNTP_dimless*Mg*param.K_Mg2NTP

    #Dimensionless concentrations for PPi
    HPPi_dimless = H*param.K_HPPi
    HMgPPi_dimless = HPPi_dimless*Mg*param.K_HMgPPi
    H2PPi_dimless = HPPi_dimless*H*param.K_H2PPi
    H2MgPPi_dimless = H2PPi_dimless*Mg*param.K_H2MgPPi
    MgPPi_dimless = Mg*param.K_MgPPi
    Mg2PPi_dimless = MgPPi_dimless*Mg*param.K_Mg2PPi

    #Dimensionless concentrations for Buffer
    HBuffer_dimless = H*10^buffer_pka
    
    #Dimensionless concentrations for Pi
    MgPi_dimless = Mg*param.K_MgPi
    HPi_dimless = H*param.K_HPi

    #Solve for Anion Concentrations
    NTP = NTPtot/(1 + HNTP_dimless + HMgNTP_dimless + MgNTP_dimless + Mg2NTP_dimless)
    PPi = PPitot/(1 + MgPPi_dimless + Mg2PPi_dimless + HPPi_dimless + HMgPPi_dimless + H2PPi_dimless + H2MgPPi_dimless)
    Buffer = Buffertot/(1+HBuffer_dimless)
    Pi = Pitot/(1+MgPi_dimless+HPi_dimless)
    OH = 10^-14/H

    #Redimensionalizing concentrations for NTP
    HNTP = HNTP_dimless*NTP
    HMgNTP = HMgNTP_dimless*NTP
    MgNTP = MgNTP_dimless*NTP
    Mg2NTP = Mg2NTP_dimless*NTP

    #Redimensionalizing concentrations for PPi
    HPPi = HPPi_dimless*PPi
    HMgPPi = HMgPPi_dimless*PPi
    H2PPi = H2PPi_dimless*PPi
    H2MgPPi = H2MgPPi_dimless*PPi
    MgPPi = MgPPi_dimless*PPi
    Mg2PPi = Mg2PPi_dimless*PPi

    #Redimensionalizing concentrations for Buffer
    HBuffer = HBuffer_dimless*Buffer
    
    #Redimensionalizing concentrations for Pi
    MgPi = MgPi_dimless*Pi
    HPi = HPi_dimless*Pi

    #Algebraic Equations
    #Change balance
    #total charge equations (correct but slower)
    negativecharge = OAc+Cl+OH+4*NTP+4*PPi+2*Pi+3*HNTP+HMgNTP+2*MgNTP+3*HPPi+HMgPPi+2*H2PPi+2*MgPPi+HPi+RNAtotalbases
    positivecharge = H+HBuffer+Na+2*Mg

    Chargebalance  = 1-(negativecharge)/(positivecharge)
    # Mg mass balance
    Mgbalance = 1 - (1/(Mgtot)) * (Mg + MgPPi + HMgPPi + MgNTP + 2*Mg2NTP + 2*Mg2PPi + HMgNTP + H2MgPPi + MgPi)

    ionspecies = (Mg,MgNTP,MgPPi,Mg2PPi)
    balances = (Chargebalance,Mgbalance)
    return (ionspecies,balances)
end

"""
    ratesmodel(param,stoich,species)

Take parameters and concentration of all species (total and free ions), return rates of processes and balances for use in DAE.

"""
function ratesmodel(param,stoich,species)
    #Species common to all reactors
    (DNAtot, RNAtot, PPitot, ATPtot, UTPtot, CTPtot, GTPtot, Mgtot, Nuctot, Pitot, H, Mg, Buffertot, RNAPtot, PPiasetot, Captot, Na, Cl, OAc) = species
    NTPtot = ATPtot + UTPtot + CTPtot + GTPtot + Captot

    #Using speciation model to calculate complex concentrations
    RNAtotalbases = RNAtot*(sum(stoich)+2)
    (ions,balances) = speciationmodel(param, RNAtotalbases, NTPtot, Mgtot, Buffertot, PPitot, Pitot, H, Mg, Na, Cl, OAc)

    (Mg,MgNTP,MgPPi,Mg2PPi) = ions

    MgNTPs = (MgNTP/NTPtot) .*[ATPtot, UTPtot, CTPtot, GTPtot] #Concentrations of MgATP, MgUTP, MgCTP, MgGTP
    
    #Calculating Transcription Rate
    if all(MgNTPs .>0)
        k_Ns = (param.k_e .*MgNTPs ./(MgNTPs .+param.K_1*(1+MgPPi/param.Ki_PPi))) .*Mg/((param.K_2+Mg))
        alpha = 1+param.k_i*sum(stoich ./ k_Ns)
        K_MD = (param.k_off+param.k_i)/param.k_on
        RNAnM = RNAPtot*1e9 # RNA in nm
        IC = 1e-9*((K_MD+RNAnM+DNAtot*alpha) - sqrt((K_MD+RNAnM+DNAtot*alpha)^2-4*RNAnM*DNAtot*alpha))/(2*alpha) #Initiation Complex concentration
        V_tr = param.k_i*IC
    else 
        V_tr = 0
    end

    #Calculating rate of nucleation and crystal growth
    S = Mg2PPi/param.Mg2PPi_eq
    if (S>1)
        B = param.B
        V_nuc = exp(-(B)/log(S)^2)
        V_precip = param.k_precip*Nuctot*log(S)
    else
        V_nuc = 0
        V_precip = 0
    end
    
    #Calculating rate of degradation of PPi by PPiase 
    if MgPPi > 0
        V_PPiase = 60*param.kPPiase*PPiasetot*MgPPi/(MgPPi+param.KMPPiase)
    else
        V_PPiase = 0
    end
    
    #Calculating rate of sequestration of DNA by nuclei
    if DNAtot > 0
        V_seq = param.k_d*(Nuctot)*DNAtot
    else
        V_seq = 0
    end

    rates = (V_seq,V_tr,V_precip,V_PPiase,V_nuc)
    return (rates,balances)
end

"""
    ivt_batch!(du,u,param,t,stoich,constantspecies)

Timestep of DAE showing time evolution of vector u by infinetesimal du.

"""
function ivt_batch!(du,u,param,t,stoich,constantspecies)
    species = vcat(u,constantspecies)
    N_all = sum(stoich)
    #Reaction Dynamics
    (rates,balances) = ratesmodel(param,stoich,species)
    (V_seq,V_tr,V_precip,V_PPiase,V_nuc) = rates
    (Chargebalance,Mgbalance) = balances
    N_A, N_U, N_C, N_G = stoich
    #Differential Equations
    #DNAtot
    du[1] =  - V_seq
    #RNAtot 
    du[2] =  V_tr
    #PPitot
    du[3] =  (N_all - 1)*(V_tr) - V_precip - V_PPiase
    #ATPtot
    du[4] =  -N_A*(V_tr)
    #UTPtot
    du[5] =  -N_U*(V_tr)
    #CTPtot
    du[6] =  -N_C*(V_tr)
    #GTPtot
    du[7] =  -N_G*(V_tr)
    #Mgtot
    du[8] =  -2*V_precip
    #Nuctot
    du[9] = V_nuc
    #Pitot
    du[10] = 2*V_PPiase

    #Algebraic Equations
    # Chargebalance
    du[11] = Chargebalance
    # Mg mass balance
    du[12] = Mgbalance
    nothing
end

"""
    runDAE_batch(params::AbstractArray{T1}, inputs; saveevery = true, stoich = SVector(231, 246, 189, 202), PPi = 1e-18, PPiase = 0.0, Cap = 0, Pi = 1e-18, Nuc = 0, RNA = 0, tol = 1e-5, init_time = 0.0) where {T1<:Number}

Take parameters and reaction inputs and run DAE model of IVT. Returns solution object.
"""
function runDAE_batch(params::AbstractArray{T1}, inputs; saveevery = true, stoich = SVector(231, 246, 189, 202), PPi = 1e-18, PPiase = 0.0, Cap = 0, Pi = 1e-18, Nuc = 0, RNA = 0, tol = 1e-6, init_time = 0.0, bufferpKa = 8.1, NaperNTP = 3.96) where {T1<:Number}
    #Generate initialization
    NTPtot = inputs.ATP+inputs.UTP+inputs.CTP+inputs.GTP
    Na = NaperNTP*NTPtot
    OAc = 2*inputs.Mg
    Cl = inputs.Buffer*((1e-8*10^bufferpKa)/(1e-8*10^bufferpKa+1)) #Cl from HCl in Buffer - remove if using HEPES   note: removed 1e-8 -1e-6

    nstatevars = 10
    nalgebraicvars = 2
    ntotalvars =  nstatevars+nalgebraicvars
    initial = zeros(T1,ntotalvars)
    initial[1:nstatevars] = [inputs.DNA,RNA,PPi,inputs.ATP,inputs.UTP,inputs.CTP,inputs.GTP,inputs.Mg,Nuc,Pi]
    NTPtot = inputs.ATP + inputs.UTP + inputs.CTP + inputs.GTP + Cap
    initial[(nstatevars+1):ntotalvars] = getfreeconcentrations(params, NTPtot, inputs.Mg, inputs.Buffer, PPi, Pi, Na, Cl, OAc)
    constantspecies = [inputs.Buffer, inputs.T7RNAP, PPiase, Cap, Na, Cl, OAc]

    #Run DAE and get solution object
    M = zeros(ntotalvars,ntotalvars)
    M[1:nstatevars,1:nstatevars] = I(nstatevars)
    f = ODEFunction(((du,u,param,t) -> ivt_batch!(du,u,param,t,stoich,constantspecies)),mass_matrix=M)
    prob_mm = ODEProblem(f,initial,(init_time,inputs.final_time),params)
    sol = solve(prob_mm,Rodas4(),abstol=tol,reltol=tol,save_everystep=saveevery)

    #If solution fails, rerun at a higher tolerance
    if (sol.retcode == ReturnCode.DtLessThanMin || sol.retcode == ReturnCode.Unstable)
        if tol>1e-12
            sol = runDAE_batch(params, inputs; saveevery = saveevery, stoich = stoich, PPi = PPi, PPiase = PPiase, Cap = Cap, Pi = Pi, Nuc = Nuc, RNA = RNA, tol = 0.1*tol, init_time = init_time)
       end
   end

    return sol
end  


