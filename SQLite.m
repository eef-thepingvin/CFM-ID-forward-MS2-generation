clear
conn=sqlite("NIST2020-HRMSMS.db","connect")
compoundTable=fetch(conn,'select * from CompoundTable');
spectrumTable=fetch(conn,'select * from SpectrumTable');

cmpIDIdx=unique(compoundTable.Name,'stable');
cmpID=1:length(cmpIDIdx);
IdxForRemoval=1:length(compoundTable.Name);

for i=1:length(compoundTable.Name)
y(i,1)=cmpID(compoundTable.Name{i} == cmpIDIdx);
end
[~,idx]=unique(y,'rows');
IdxForRemoval(idx)=[]; %Indices for row removal in compoundTable
x=unique(y);

spectrumIDs=array2table(int64(y))
spectrumIDs.Properties.VariableNames={'NewIDs'}
compoundIdx=array2table(int64(IdxForRemoval'))
compoundIdx.Properties.VariableNames={'removeIdx'}
newCompoundIdx=array2table(int64(x))
newCompoundIdx.Properties.VariableNames={'NewIdx'}

sqlwrite(conn,'CompoundRowRemove',compoundIdx)
sqlwrite(conn,'newIDX',newCompoundIdx)
sqlwrite(conn,'SpectrumIDs',spectrumIDs)