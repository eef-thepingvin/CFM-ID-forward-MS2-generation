using DataFrames
using CSV
using PyCall
using Conda

#Definitions
filename="cfm-id-output-neg.msp";
#mode="positive"
mode="negative"

##Installing rdkit 
#Conda.add("rdkit"; channel="conda-forge") #requires conda v. 23.1.0
##If errors occurs, i.e. wrong version, use following commands in conda
#conda activate ~/.julia/conda/3/x86_64/
#conda install conda=23.1.0

## Installing cfm-id python package - DOES NOT CURRENTLY WORK
#Conda.add("llvmlite";channel="numba")
#Conda.add("numba")
#Conda.pip("install","numba")
#Conda.pip("install","matchms")
#Conda.pip("install","cfm-id") #Might fail at Numba installation due to missing Microsoft Visual C++ ver. 14.0 or greater. Install this first. #Might need to downgrade Python to v.3.10: conda install python=3.10
#Conda.pip_interop(true)
#Conda.pip("install","git+git://github.com/hcji/PyCFMID@master")

#Loads RDKit package(s)
rdk=pyimport("rdkit.Chem") #Loads RDKit module from Python/Conda
rdkRDLog=pyimport("rdkit.RDLogger")
rdkDes=pyimport("rdkit.Chem.Descriptors") #Loads Descriptors module from RDKit to calculate monoisotopic mass from SMILES
rdMolDesc=pyimport("rdkit.Chem.rdMolDescriptors")  #Loads rdMolDescriptors module from RDKit to calculate molecular formula from SMILES
rdkInchi=pyimport("rdkit.Chem.inchi") #Loads inchi module from RDKit to calculate InChI and InChIKeys from SMILES
rdkRDLog.DisableLog("rdApp.*")

dataList=readlines(filename)
println("Step 1 of 5: Performing string replacement of $(length(dataList)) lines")

#String replacement
dataList=replace.(dataList,"_" => " ");
dataList=replace.(dataList,"+ve in-silico MS/MS by CFM-ID 4.4.7 for " => "");
dataList=replace.(dataList,"-ve in-silico MS/MS by CFM-ID 4.4.7 for " => ""); 
dataList=replace.(dataList,"Comment:" => "Collision_energy:");
dataList=replace.(dataList,"Energy0" => "10V\nIonMode:\nFormula: \nAdduct:\nPrecursorMz: \nInChIKey: \nInChI: \nComments:");
dataList=replace.(dataList,"Energy1" => "20V\nIonMode:\nFormula: \nAdduct:\nPrecursorMz: \nInChIKey: \nInChI: \nComments:");
dataList=replace.(dataList,"Energy2" => "40V\nIonMode:\nFormula: \nAdduct:\nPrecursorMz: \nInChIKey: \nInChI: \nComments:");
dataList=replace.(dataList,"Smiles/Inchi:" => "SMILES: ");
dataList=replace.(dataList,"Comments:" => "Comments: In-silico MS/MS prediction by CFM-ID 4.4.7");
    if contains(mode,"positive")
        dataList=replace.(dataList,"IonMode:" => "IonMode: positive");
        dataList=replace.(dataList,"Adduct:" => "Precursor_type: [M+H]+"); 
    else
        dataList=replace.(dataList,"IonMode:" => "IonMode: negative");
        dataList=replace.(dataList,"Adduct:" => "Precursor_type: [M-H]-"); 
    end
dataList=dataList.*"\n"

#Saving and reimporting list of lines
write("temp.msp",join(dataList)) #generates temporary data file;
dataList=readlines("temp.msp");
rm("temp.msp") #deletes temporary data file

#Extracting SMILES
smilesIdx=findall( x -> occursin("SMILES: ", x), dataList)
smiles=dataList[smilesIdx]
smiles=replace.(smiles,"SMILES: " => "")

#Calculates monoisotopic mass
println("Step 2 of 5: Calculating $(length(smiles)) molblocks")
massVector=ones(length(smiles))
identifierMatrix=string.(ones(length(smiles),3))
for i in eachindex(smiles)
    massVector[i,1]=rdkDes.ExactMolWt(rdk.MolFromSmiles(smiles[i])) #Calculates monoisotopic mass
    identifierMatrix[i,1]=rdMolDesc.CalcMolFormula(rdk.MolFromSmiles(smiles[i])) #Calculates formula
    identifierMatrix[i,2]=rdkInchi.MolToInchiKey(rdk.MolFromSmiles(smiles[i])) #Calculates InChiKey
    identifierMatrix[i,3]=rdkInchi.MolToInchi(rdk.MolFromSmiles(smiles[i])) #Calculates InChiKey
    println("Step 2 of 5: Calculating $i of $(length(smiles)) molblocks")
end

#Adds/subtracts the mass of a proton and converts to string
protonMass=1.00727647;
if contains(mode,"negative")==1
    protonMass=-protonMass;
else
end
println("Step 3 of 5: Calculating precursor masses")
massVector=massVector .+ protonMass
massVector=string.(round.(massVector,digits=5))

#Inserts calculated monoisotopic masses, formula, and InChIKey
println("Step 4 of 5: Inserting $(length(smiles)) identifiers")
massIdx=findall( x -> occursin("PrecursorMz", x), dataList)
formulaIdx=findall( x -> occursin("Formula", x), dataList)
inchikeyIdx=findall( x -> occursin("InChIKey", x), dataList)
inchiIdx=findall( x -> occursin("InChI:", x), dataList)
for i in eachindex(massIdx)
    dataList[massIdx[i]] = dataList[massIdx[i]]*massVector[i]
    dataList[formulaIdx[i]] = dataList[formulaIdx[i]]*identifierMatrix[i,1]
    dataList[inchikeyIdx[i]] = dataList[inchikeyIdx[i]]*identifierMatrix[i,2]
    dataList[inchiIdx[i]] = dataList[inchiIdx[i]]*identifierMatrix[i,3]
end

# Find errors between assignmed peak numbers and actual peak numbers from errorenous CFM-ID computation and replaces with correct number of peaks
println("Step 5 of 5: Replacing $(length(smiles)) feature numbers")
a=findall(x -> occursin("Num peaks:", x), dataList) #Index all lines defining number of features in each peak
#aa=parse.(Int64,filter.(isdigit,dataList[a])) #Extracts number of defined features in each peak and converts to Int64
b=findall(isempty.(dataList)) .- 1 #Index each empty line following features
c=string.(b.-a) #Finds true number of peaks and converts to string
for i in eachindex(a)
    dataList[a[i]]="Num peaks: "*c[i] #Replaces peak number values with correct values
end

#Prepares datalist for writelines in Julia by join-function
dataList=dataList.*"\n"

write(replace(filename,".msp" => "_Julia_cleanup_DONE.msp"),join(dataList));
println("Convertion done!")