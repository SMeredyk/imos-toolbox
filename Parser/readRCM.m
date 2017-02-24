function sample_data = readRCM( filename, mode )
% readRCM reads a data file retrieved from an RCM4,7,11 data logger.
%
% This function is able to read in a single text file (tab delimited) retrieved from an 
% RCM4,7,11 data logger exported by the 5059 Data Reading software. 
% The pressure data is returned in a sample_data struct. 
% 
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Header Example:
% Model: RCM4 / RCM7 / RCM11
% ProductNumber:	4430 (possibly blank)
% SerialNumber:	30
%
% Author : 		 Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributors : Pascal_Guillot@uqar.ca (UQAR - Canada), Guillaume Galibert <guillaume.galibert@utas.edu.au>			
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

% open the file, and read in the data
try
    fid = fopen(filename, 'rt');
    rawText = textscan(fid, '%s', 'Delimiter', ''); % comma delimited text file
    fclose(fid);
    
    [header, iData] = readHeader(rawText{1});
    data = readData(filename, iData);
catch e
    if fid ~= -1, fclose(fid); end
    rethrow(e);
end % end of File read

% copy all of the information over to the sample data struct
sample_data = struct;

sample_data.toolbox_input_file              = filename;
sample_data.meta.instrument_make            = 'Aanderaa';
sample_data.meta.instrument_model           = header.Model;
sample_data.meta.instrument_serial_no       = header.SerialNumber;
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
data = rmfield(data, {'TIME','REF'});		
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

function [header, iData] = readHeader(rawText)
%READHEADER Reads the header , including channel coefficients from the file.

  header 	= struct;
  iData 	= [];
  
  startHeader 	= 'Model:RCM11';      
  endHeader 	= '[Data]';   
  fmtHeader  	= '%s%s';
  delimHeader 	= ':';
  
  
  iStartHeader = find(strcmp(startHeader, rawText));
  iEndHeader = find(strcmp(endHeader, rawText))-1;
  iData = iEndHeader + 1;
  
  headerCell 	= rawText(iStartHeader:iEndHeader);
  nFields 		= length(headerCell);
    
  for i=1:nFields
      tuple = textscan(headerCell{i}, fmtHeader, 'Delimiter', delimHeader);
      header.(tuple{1}{1}) = tuple{2}{1};
  end
end % end of readHeader function


function data = readData(filename, iData)
%READDATA Reads the sample data from the file.

  data = struct;
  dataDelim = ',';	% comma delimited data
  
  fid = fopen(filename, 'rt');
  params = textscan(fid, '%s', 1, 'HeaderLines', iData, 'Delimiter', '');   % iData passed the header position to readData
  params = params{1};
  iParams = strfind(params, ',');
  nParams = length(iParams{1});
  paramsFmt = repmat('%s', 1, nParams);
  params = textscan(params{1}, paramsFmt, 'Delimiter', dataDelim);
  dataFmt = ['%s', repmat('%f', 1, nParams-1)];
  values = textscan(fid, dataFmt, 'Delimiter', dataDelim);
  fclose(fid);
  
  for i=1:nParams
      switch params{i}{1}
                  
                  %Date Time (mm.dd.yyyy hh:mm:ss) 
                  case 'Time tag (Gmt)', 
                    name = 'TIME';
                    data.TIME.values = datenum(values{i}, 'mm.dd.yyyy HH:MM:SS');
					data.TIME.comment = ['Date Time (mm.dd.yyyy hh:mm:ss)'];
					
                  %Reference Parameter (unitless) - not used by Toolbox
                  case 'Reference', 
                    name = 'REF';
                    data.REF.comment = [''];
				
				  %Turbidity (Ftu) - NTU and FTU have similar values ,
				  %though calibration method is different chemicals
				  %(Forazine , Nephelometric)
                  
				  case 'Turbidity(Ftu)', 
				    name = 'TURBF';
                    data.TURBF.values = values{i};
					data.TURBF.comment = ['Turbidity data '...
					  'computed from bio-optical sensor raw counts measurements. The '...
					  'turbidity sensor is equipped with a 880nm peak wavelength LED to irradiate and a '...
					  'photodetector paired with an optical filter which measures everything '...
					  'that backscatters in the region of 650nm to 1000nm.']';

                  case 'Turbidity(NTU)', 
				    name = 'TURB';
                    data.TURB.values = values{i};
					data.TURBF.comment = ['Turbidity data '...
					  'computed from bio-optical sensor raw counts measurements. The '...
					  'turbidity sensor is equipped with a 880nm peak wavelength LED to irradiate and a '...
					  'photodetector paired with an optical filter which measures everything '...
					  'that backscatters in the region of 650nm to 1000nm.']';
                  
                  % Two pressure resolution options MPa or kPa, depending on .cdb file
                  % export from 5059 software
                  
				  %Pressure (MPa) = 100-1*(dBarr)
                  case 'Pressure(MPa)', 
                     name = 'PRES';
                     data.PRES.values = (values{i}*100);   
					 data.PRES.comment = ['Pressure data converted from MPa to dBarr for toolbox'];	
				
                  %Pressure (kPa) = 10-1*(dBarr) 
                  case 'Pressure(kPa)', 
                     name = 'PRES';
                     data.PRES.values = (values{i}*10);   
					 data.PRES.comment = ['Pressure data converted from kPa to dBarr for toolbox'];	

				  %Temperature (Celsius degree)
                  case 'Temperature(DegC)', 
                    name = 'TEMP';
                    data.TEMP.values = values{i};
					data.TEMP.comment = ['Each sensor on the data logger comes with its own temperature ref,'....
										' though only temperature sensor is used'];
                                    
				  %Conductivity (mS/cm) = 10-1*(S/m)
                  case 'Conductivity(mS/cm)',
                      name = 'CNDC';
                      data.CNDC.values = (values{i})/10; 
                      data.CNDC.comment = ['converted from mS/cm to S/m by Toolbox'];
					  
				  %DO (uM/L)
                  case 'O2Concentration(uM)', 
                    name = 'DOX1';
                    data.DOX1.values = values{i};
                    data.DOX1.comment = ['Dissolved Oxygen Concentration (uM)per Litre'];
                        					  
                  % it's better to let the toolbox calculate the Depth (gsw) 
                  %Depth (m)	Latitude corrected in Aanderaa software
                  % case 'Depth(m)', 
					% name = 'DEPTH';
                    %data.DEPTH.values = values{i};
					% data.DEPTH.comment = ['Derived SeaWater Depth Calculated by Aanderaa Seaguard Studio'...
						%					' and latitude corrected'];
						
                  %Density(kg/m^3)	 -  its better to let the toolbox derive this
                   % case 'Density(kg/m^3)', 
					%	name = 'DENS';
                    %   data.DENS.values = values{i};
					%	data.DENS.comment = ['Derived SeaWater Depth Calculated by Aanderaa Seaguard Studio'...
					%						' and latitude corrected'];  
                  
                 %% Single-Point Current Meter - Very Basic current Measurements - assumed no mag. decl. applied
                 
                 %Absolute Speed (cm/s) = 100-1*(m/s)
                  case 'AbsSpeed(cm/s)', 
                    name = 'CSPD';
                    data.CSPD.values = (values{i})/100; 
					data.CSPD.comment = ['Absolute SeaWater Velocity m/s'];
					
				%Absolute Water Direction (Deg.M)
                  case 'Direction(Deg.M)', 
                    name = 'CDIR';
                    data.CDIR.values = values{i};
					data.CDIR.comment = ['Direction of the SeaWater Velocity'];
				
      end % end of Switch
  end % end of For loop
end % end of readData function
