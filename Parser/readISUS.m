function sample_data = readISUS(filename, mode )
%readISUS Parses several data files retrieved from a Satlantic nitrate sensor ISUS V2 /V3 
%
% This function is able to read a concatenated dataset of .dat files retrieved from an ISUS unit.
% The data is recovered via the ISUSCom 2.1.6 software and converted into a
% single CSV file via SatCon 1.5.5.2. The Date and Time
% columns have been concatenated and formatted in DATETIME as dd-mm-yyyy hh:mm:ss.
%
% Inputs:
%   filename  	- Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 		 Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributors : Guillaume Galibert <guillaume.galibert@utas.edu.au>			
%
% Copyright (c) 2017, Amundsen Science & ArcticNet
% http://www.amundsen.ulaval.ca/
% http://www.arcticnet.ulaval.ca/
% All rights reserved.
%
% Copyright (c) 2016, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS).
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%     * Redistributions of source code must retain the above copyright notice, 
%       this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in the 
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the AODN/IMOS nor the names of its contributors 
%       may be used to endorse or promote products derived from this software 
%       without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%
% ensure that there is one or two arguments
narginchk(1,2);

if ~ischar(filename), error('filename must contain a string'); end

[~, name, ~] = fileparts(filename);
%extracting filename parts 
unitInfo = textscan(name, '%s %s %d', 'Delimiter', '_');
%unitMake = char(unitInfo{1});
%unitModel = char(unitInfo{2});

data = readData(filename);   
 
% copy all of the information over to the sample data struct
sample_data = struct;

sample_data.toolbox_input_file              = filename;
sample_data.meta.instrument_make            = char(unitInfo{1}); 
sample_data.meta.instrument_model           = char(unitInfo{2});
sample_data.meta.instrument_serial_no       = unitInfo{3}; 
sample_data.meta.instrument_sample_interval = median(diff(data.TIME.values*24*3600)); 
sample_data.meta.featureType                = mode;

sample_data.dimensions = {};
sample_data.variables  = {};

sample_data.dimensions{1}.name              = 'TIME';
sample_data.dimensions{1}.typeCastFunc      = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{1}.name, 'type')));
sample_data.dimensions{1}.data              = sample_data.dimensions{1}.typeCastFunc(data.TIME.values);

sample_data.variables{end+1}.name           = 'TIMESERIES';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(1);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'LATITUDE';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'LONGITUDE';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];
sample_data.variables{end+1}.name           = 'NOMINAL_DEPTH';
sample_data.variables{end}.typeCastFunc     = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
sample_data.variables{end}.data             = sample_data.variables{end}.typeCastFunc(NaN);
sample_data.variables{end}.dimensions       = [];

% copy variable data over
data = rmfield(data, 'TIME');
fields = fieldnames(data);

for k = 1:length(fields)
    name = fields{k};
    
    % dimensions definition must stay in this order : T, Z, Y, X, others;
    % to be CF compliant
    sample_data.variables{end+1}.dimensions = 1;
    sample_data.variables{end}.name         = name;
    sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
    sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc(data.(fields{k}).values);
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
    sample_data.variables{end}.comment      = data.(fields{k}).comment;
end
end %end of main function


function data = readData(filename)
%READDATA Reads the sample data from the file.

  data = struct;
  dataDelim = ',';	% comma delimited data 
  
  fid = fopen(filename, 'rt');
  params = textscan(fid, '%s', 1, 'Delimiter', '');   
  params = params{1};
  iParams = strfind(params, ',');
  nParams = length(iParams{1})+1; %
  paramsFmt = repmat('%s', 1, nParams);
  params = textscan(params{1}, paramsFmt, 'Delimiter', dataDelim); 
  dataFmt = ['%s', repmat('%f', 1, nParams-1)];  
  values = textscan(fid, dataFmt, 'Delimiter', dataDelim); 
  fclose(fid);
  
% other variables exported by SatCon and avaiable for analysis
%   T_SPEC	T_LAMP	HUMIDITY    VOLT_12	VOLT_5	REF_AVG	REF_STD
%   SW_DARK	SPEC_AVG	UV(189.17) - all wavelengths

  for i=1:nParams
      switch params{i}{1}  
				  				
				%DateTime (dd-mm-yyyy HH:MM:SS)
                  case 'DATETIME',
                    name = 'TIME';
                    data.TIME.values = datenum(values{i},'dd-mm-yyyy HH:MM:SS');			
                    data.TIME.comment = ['dd-mm-yyyy HH:MM:SS'];  				
										  
				  case 'NITRATE', 
				    name = 'NTRA';
                    data.NTRA.values = values{i};
					data.NTRA.comment = ['mole_concentration_of_nitrate_in_sea_water : [mole l-1]']';
				
				case 'T_INT', 
				    name = 'DRYT';
                    data.DRYT.values = values{i};
					data.DRYT.comment = ['degrees celcius - dry internal temp of nitrate sensor lamp']';
                    
                case 'T_LAMP', 
				    name = 'TEMP';
                    data.TEMP.values = values{i};
					data.TEMP.comment = ['degrees celcius - dry internal temp of nitrate sensor']';
                
				case 'VOLT_MAIN', 
				    name = 'VOLT';
                    data.VOLT.values = values{i};
					data.VOLT.comment = ['Voltage main of nitrate sensor for diagnostic purposes']';  
                    
      end % end of Switch
  end % end of For loop
end % end of readData function