"""
    importakamadata()

Process akama data files for use in model fitting."""
function importakamadata()
    #Importing Akama Intial Rate Data
    df = CSV.read("data/akama_data/akama_initial_rate.csv", DataFrame)
    rawdata = Matrix(df)
    #Convert units of Akama rate data
    processedratedata = zeros(size(rawdata)[1],4) 

    # Columns of data matrix: NTP (tot), Mg(tot), Mg (free), MgNTP, RNA, RNA St dev

    processedratedata[:,1] = 4 * 0.001 .* rawdata[:,1]     # NTP mM (each) to M (total)
    processedratedata[:,2] = 0.001 .* rawdata[:,2]         # Mg mM to M

    processedratedata[:,3] = 5*1e-6 .* rawdata[:,3]        # RNA yield from uM/min (over 5 min) to M
    processedratedata[:,4] = 5*1e-6 .* rawdata[:,4]        # RNA yield from uM/min (over 5 min) to M
    filteredratedata = processedratedata[(processedratedata[:,1] .!= 0),:]#Remove [NTP] = 0 points
    initialrateconcentrationinputs = filteredratedata[:,1:2]
    initialrateyieldoutputs = filteredratedata[:,3]
    initialrateyieldstdev = filteredratedata[:,4];

    #Importing Akama Mg titration
    df = CSV.read("data/akama_data/akama_mg2ppi_solubility.csv", DataFrame)
    rawdatamat = Matrix(df)
    #Convert units of Akama rate data
    processeddata = similar(rawdatamat)
    processeddata[:,1] = 4 * 0.001 .* rawdatamat[:,1]      #NTP from mM (each) to M (total)
    processeddata[:,2] = 0.001 .* rawdatamat[:,2]          # PPi from mM to M
    processeddata[:,3] = 0.001 .* rawdatamat[:,3]          # Mg from mM to M
    processeddata[:,4] = 0.001 .* rawdatamat[:,4]          # Mg st dev  from mM to M

    Mgconcentrationinputs = processeddata[:,1:2]
    Mgoutputs = processeddata[:,3]
    Mgstdev = processeddata[:,4];

    #Importing Akama Rate Data
    df = CSV.read("data/akama_data/akama_reaction_trajectories.csv", DataFrame)
    rawdatamat = Matrix(df)
    #Convert units of Akama rate data
    processeddata = similar(rawdatamat)
    processeddata[:,1] = (1/60) .* rawdatamat[:,1]         # Time from min to hr 
    processeddata[:,2] = 1e-6 .* rawdatamat[:,2]           # T7RNAP from uM to M
    processeddata[:,3] = 4 * 0.001 .* rawdatamat[:,3]      # NTP from mM (each) to M (total)
    processeddata[:,4] = 0.001 .* rawdatamat[:,4]          # Mg from mM to M
    processeddata[:,5] = 1e-6 .* rawdatamat[:,5]           # RNA from uM to M
    processeddata[:,6] = 1e-6 .* rawdatamat[:,6]           # RNA from std from uM to M
    processeddata[:,7] = 0.001 .* rawdatamat[:,7]          # PPi from mM to M
    processeddata[:,8] = 0.001 .* rawdatamat[:,8]          # PPi st dev from mM to M

    filteredratedata = processeddata[(processeddata[:,1] .!= 0),:];#Remove t = 0 points - trivially fit all models
    concentrationinputs, timeinputs, RNAyieldoutputs, PPiyieldoutputs, RNAstdev, PPistdev = transformakamadata(filteredratedata)

    
    return (; concentrationinputs, timeinputs, RNAyieldoutputs, PPiyieldoutputs, RNAstdev, PPistdev, Mgconcentrationinputs, Mgoutputs, Mgstdev, initialrateconcentrationinputs, initialrateyieldoutputs, initialrateyieldstdev)
end

"""
    transformakamadata(akamadata)

Transform imported CSV into format for data structure."""
function transformakamadata(akamadata)
    numberofindependentconditions = Int64(size(akamadata)[1]/9)
    concentrationinputs = zeros(numberofindependentconditions,3)
    timeinputs = zeros(numberofindependentconditions,9)
    RNAyieldoutputs = zeros(numberofindependentconditions,9)
    PPiyieldoutputs = zeros(numberofindependentconditions,9)
    RNAstdev = zeros(numberofindependentconditions,9)
    PPistdev = zeros(numberofindependentconditions,9)
    RNAstdevmin = sqrt((norm(akamadata[:,6])^2)/size(akamadata)[1])
    PPistdevmin = sqrt((norm(akamadata[:,8])^2)/size(akamadata)[1])
    for i in 1:Int64(size(akamadata)[1]/9)
        concentrationinputs[i,:] = akamadata[i*9-8,2:4]
        timeinputs[i,:] = akamadata[i*9-8:i*9,1]
        RNAyieldoutputs[i,:] = akamadata[i*9-8:i*9,5]
        PPiyieldoutputs[i,:] = akamadata[i*9-8:i*9,7]
        RNAstdev[i,:] = max.(akamadata[i*9-8:i*9,6],RNAstdevmin)
        PPistdev[i,:] = max.(akamadata[i*9-8:i*9,8],PPistdevmin)
    end
    return concentrationinputs, timeinputs, RNAyieldoutputs, PPiyieldoutputs, RNAstdev, PPistdev
end

"""
    importakamaPPiasedata()

Import akama data for PPiase use (main text figure 4)."""
function importakamaPPiasedata()
    #Importing Akama Rate Data
    df = CSV.read("data/akama_data/akama_PPiase.csv", DataFrame)
    rawdatamatPPiase = Matrix(df)
    #Convert units of Akama rate data
    processedPPiasedata = similar(rawdatamatPPiase)
    processedPPiasedata[:,1] = (1/60) .* rawdatamatPPiase[:,1]         # Time from min to hr 
    processedPPiasedata[:,2] = 1e-6 .* rawdatamatPPiase[:,2]           # RNAfrom uM to M
    processedPPiasedata[:,3] = 1e-6 .* rawdatamatPPiase[:,3]           # RNAfrom uM to M
    processedPPiasedata[:,4] = 0.001 .* rawdatamatPPiase[:,4]          # PPi from mM  to M 
    processedPPiasedata[:,5] = 0.001 .* rawdatamatPPiase[:,5]          # PPi from mM  to M 
    processedPPiasedata[:,6] = 1e-6 .* rawdatamatPPiase[:,6]           # RNAfrom uM to M
    processedPPiasedata[:,7] = 1e-6 .* rawdatamatPPiase[:,7]           # RNAfrom uM to M
    processedPPiasedata[:,8] = 0.001 .* rawdatamatPPiase[:,8]          # PPi from mM  to M 
    processedPPiasedata[:,9] = 0.001 .* rawdatamatPPiase[:,9]
    processedPPiasedata
end