function y = IpsExtractCsvV650VertCat (filename)
%% This function reads all CSV files and vertically 
% concatenates them into a single table of all ice draft and ice draft
% error values for the IMOS process chain.




ds = datastore(pwd);
mydata = ds.readall;

end