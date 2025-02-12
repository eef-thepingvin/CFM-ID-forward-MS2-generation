using PyCall
using Conda
using DataFrames
using CSV

#Conda.add("rdkit"; channel="conda-forge") #requires conda v. 23.1.0

#If errors occure, i.e. wrong version, use following commands in conda
#conda activate ~/.julia/conda/3/x86_64/
#conda install conda=23.1.0

rdk=pyimport("rdkit.Chem")
rdkSalt=pyimport("rdkit.Chem.SaltRemover")
rdkDes=pyimport("rdkit.Chem.Descriptors")
rdkInchi=pyimport("rdkit.Chem.inchi")
rdMolDesc=pyimport("rdkit.Chem.rdMolDescriptors")
rdMolStandardize=pyimport("rdkit.Chem.MolStandardize.rdMolStandardize")
rdkmf=pyimport("rdkit.Chem.rdmolfiles")
rdkRDLog=pyimport("rdkit.RDLogger")
rdkRDLog.DisableLog("rdApp.*")

#List of SMILES
importSMILES=CSV.read("smilesInput.csv",DataFrame)
SMILES=String.(importSMILES.SMILES)
itr=length(SMILES)
#numericMolblock=DataFrame(Mass = ones(itr))
stringMolblock=[SMILES SMILES SMILES SMILES SMILES ones(itr) ones(itr) ones(itr)]

#Salt Definitions
remover=rdkSalt.SaltRemover()#(defnData="[Na,Fe,Ca]")#(defnData="CC")

#SMILES cleanup of defined salts
for i in eachindex(SMILES)
try
        #stringMolblock[i,2]=rdk.MolToSmiles(remover.StripMol(rdk.MolFromSmiles(SMILES[i]),dontRemoveEverything=false,sanitize=true))
        stringMolblock[i,2]=rdk.MolToSmiles(rdMolStandardize.ChargeParent(rdk.MolFromSmiles(SMILES[i]))) #Removes salt and neutralizes SMILES from 
        stringMolblock[i,3]=rdMolDesc.CalcMolFormula(rdk.MolFromSmiles(stringMolblock[i,2])) #Calculates formula
        stringMolblock[i,4]=rdkInchi.MolToInchi(rdk.MolFromSmiles(stringMolblock[i,2])) #Calculates InChi
        stringMolblock[i,5]=rdkInchi.MolToInchiKey(rdk.MolFromSmiles(stringMolblock[i,2])) #Calculates InChiKey
        stringMolblock[i,6]=rdk.GetFormalCharge(rdk.MolFromSmiles(stringMolblock[i,2])) #Calculates formal charge of molecule from neutralized SMILES
        stringMolblock[i,7]=round(rdkDes.ExactMolWt(rdk.MolFromSmiles(stringMolblock[i,2])),digits=8) #Calculates monoisotopic mass
        stringMolblock[i,8]=rdk.Mol.GetNumAtoms(rdk.AddHs(rdkmf.MolFromSmiles(stringMolblock[i,2]))) #Calculates number of atoms in structure (important for CFM-ID as n_atoms < 251        println("Calculating $i of $itr molblocks")
        println("Calculating $i of $itr molblocks")
catch e
        continue
end
end

#smiles = "CCC.CC"
#fragment = max(smiles.split('.'), key=len)
#print (fragment)

#####THIS PART FOR CDv3.2 or CDv3.3 masslists ONLY
#for i=1:itr
#    stringMolblock[i,5]=replace(stringMolblock[i,4],"\n" => ";") #this is required for CD masslists - as they for some reason messed up the format since CD v3.2
#end

cd("C:/Users/emefro/Desktop/Files/PhD/Julia/RDKit")
CSV.write("output-saltremover.csv",DataFrame(SMILES_input = stringMolblock[:,1], SMILES_neutralised=stringMolblock[:,2], Formula=stringMolblock[:,3], InChI=stringMolblock[:,4], InChIKey=stringMolblock[:,5], Charge=stringMolblock[:,6], MonoisotopicMass=stringMolblock[:,7], NumberOfAtoms=stringMolblock[:,8]))
