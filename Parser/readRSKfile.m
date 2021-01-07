function sample_data = readRSKfile(fname, mode)
%% Function Description
% readRSKfile.m Parses a data file retrieved from an RBR logger 
% that was exported from the new (fall 2020) Ruskin(c) v2.12 software.
% this script will extract the data needed to process an .rsk file
% without having to manually export a Ruskin Zip folder.
% It uses existing open-source Matlab scripts from RBR to speed-up coding
% though it requires the RSKTools folder to be added as a path in Matlab.
% The RSKTools directory is within the Utils directory of the IMOS toolbox.
%
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('timeSeries').
%                 *profiling mode is not supported at this time (jan, 2021)
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 		Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributor : Guillaume Galibert <guillaume.galibert@utas.edu.au>,
%               and selected RBR RSKTools scripts (2020).
%				
%
% Copyright (c) 2021, Amundsen Science & ArcticNet
% http://www.amundsen.ulaval.ca/
% http://www.arcticnet.ulaval.ca/
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
% End of Function Description

  narginchk(2,2);
  
  if ~ischar(fname)  
    error('filename must be a string'); 
  end
  
  % open the file, and read in the header and data
  try 
    
    %% opening an .rsk file using RBR RSKtools scripts
    rsk = RSKopen(fname); 
    rsk = RSKreaddata(rsk);
    RSK = RSK2MAT(rsk);

  catch e
    if fname ~= -1, fclose(fname); end
    rethrow(e);
  end
  
  %% importing extracted data and metadata into the sample_data struct
  % to be IMOS compliant
  data = struct;
  % adding time to data struct
  data.TIME.values = rsk.data.tstamp;
  
  % copy all of the information over to the sample data struct
  sample_data = struct;

    %extracting name parts from RSK.name
    unitInfo = textscan(RSK.name, '%s %s %s %d', 'Delimiter', ' ');
  
    sample_data.toolbox_input_file                = fname;
    sample_data.meta.instrument_make              = char(unitInfo{1});
    sample_data.meta.instrument_model             = char(unitInfo{2});
    sample_data.meta.instrument_firmware          = char(unitInfo{3});
    sample_data.meta.instrument_serial_no         = double(unitInfo{4});
    sample_data.meta.instrument_sample_interval   = rsk.average.repetitionPeriod ./ 60000; %interval in minutes
    sample_data.meta.featureType                  = mode; %assumed timeSeries

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

 
 %% transfering the RSK struct into the data struct
 for i=1:length(RSK.channelnames)
      switch RSK.channelnames{i,1}
                  		
				  %Turbidity (Ftu) - NTU and FTU have similar values ,
				  %though calibration method is different chemicals
				  %(Forazine , Nephelometric)
                  
				  case 'Turbidity'
				    name{i} = 'TURB';
                    data.TURB.values = RSK.data(:,i); 
					data.TURB.comment = ['NTU'];
			
                  %Pressure(dBarr) 
                  case 'Pressure' 
                     name{i} = 'PRES';
                     data.PRES.values = RSK.data(:,i);   
					 data.PRES.comment = ['dbar'];	

				  %Temperature (Celsius degree)
                  case 'Temperature' 
                    name{i} = 'TEMP';
                    data.TEMP.values = RSK.data(:,i); 
					data.TEMP.comment = ['Degrees Celius'];
                                    
				  %Conductivity (mS/cm) = 10-1*(S/m)
                  case 'Conductivity'
                      name{i} = 'CNDC';
                      data.CNDC.values = (RSK.data(:,i))/10; 
                      data.CNDC.comment = ['converted from mS/cm to S/m for Toolbox'];
					  
				  %DO (uM/L)
                  case 'O2Concentration' 
                    name{i} = 'DOX1';
                    data.DOX1.values = RSK.data(:,i);
                    data.DOX1.comment = ['Dissolved Oxygen Concentration (uM)per Litre'];           
				
      end % end of Switch
  end % end of For loop
             
%%%
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

%house cleaning
clear rsk RSK fname data name fields k i
  
end % end of script / function  