function sample_data = readRDCP(filename, mode)
% readRDCP reads a data file retrieved from a Seaguard RDCP data logger.
%
% This function is able to read in a single text file (tab delimited) retrieved from an 
% Seaguard data logger exported by the Seaguard RDCP software. 
% The pressure data is returned in a sample_data struct. 
% 
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Model::RDCP
% SerialNumber::270
% RefReading::593
% BlankDist::1
% CellSize::4
% NumCells::30
% CellOverlap::50
% DistToFirstBin::3
% RecordInterval::3600
% StartTime::14.08.2007 01:00
% 
% [Data]
%
% Author :      Guillaume Galibert <guillaume.galibert@utas.edu.au>
% Contributors: Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
%               Pascal_Guillot@uqar.ca (UQAR - Canada)			
%
% Copyright (c) 2010, eMarine Information Infrastructure (eMII) and Integrated 
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
%     * Neither the name of the eMII/IMOS nor the names of its contributors 
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
    rawText = textscan(fid, '%s', 'Delimiter', ''); % tab delimited text file
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
sample_data.meta.head                       = 600;
sample_data.meta.instrument_sample_interval = median(diff(data.TIME.values*24*3600));
sample_data.meta.beam_angle                 = 25;  
sample_data.meta.binSize                    = header.CellSize; 
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

function [header, iData] = readHeader(rawText)
%READHEADER Reads the header , including channel coefficients from the file.

  header 	= struct;
  iData 	= [];
  
  startHeader 	= 'Model:RDCP';	
  endHeader 	= '[Data]';    
  fmtHeader  	= '%s%s';
  delimHeader 	= '::';
  
  
  iStartHeader 	= find(strcmp(startHeader, rawText));
  iEndHeader 	= find(strcmp(endHeader, rawText));
  iData 		= iEndHeader + 1;   
  
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
  dataDelim = ',';	% tab delimited data...or could be space.... to be tested
  
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
                  %Date Time (dd.mm.yy HH:MM:SS)
                  case 'Time tag (Gmt)'
                    data.TIME.values = datenum(values{i}, 'dd.mm.yyyy HH:MM:SS');
                    data.TIME.comment = ['dd.mm.yyyy HH:MM:SS'];		
              
                  %Battery (Volts)
                  case 'BatteryVoltage', 
                    name = 'VOLT';
					data.VOLT.comment = ['Volts'];
				
				  %Turbidity (Ftu) - NTU and FTU have similar values ,
				  %though calibration method is different chemicals
				  %(Forazine , Nephelometric)
				  
% 				  case 'Turbidity(Ftu)', 
% 				    name = 'TURBF';
% 					data.TURBF.comment = ['Turbidity data '...
% 					  'computed from bio-optical sensor raw counts measurements. The '...
% 					  'turbidity sensor is equipped with a 880nm peak wavelength LED to irradiate and a '...
% 					  'photodetector paired with an optical filter which measures everything '...
% 					  'that backscatters in the region of 650nm to 1000nm.']';
% 
%                   case 'Turbidity(NTU)', 
% 				    name = 'TURB';
% 					data.TURB.comment = ['Turbidity data '...
% 					  'computed from bio-optical sensor raw counts measurements. The '...
% 					  'turbidity sensor is equipped with a 880nm peak wavelength LED to irradiate and a '...
% 					  'photodetector paired with an optical filter which measures everything '...
% 					  'that backscatters in the region of 650nm to 1000nm.']';
%                   
				  % Two pressure options MPa or kPa, depending on .cdb file
                  % export from 5059 software
                  
				  %Pressure (MPa) = 100-1*(dBarr)
                  %case 'Pressure(MPa)', 
                  %   name = 'PRES';
                  %   data.(fields{k}) = data.(fields{k})/100;  
					% data.PRES.comment = ['Pressure data in kPa converted to dBarr for toolbox'];	
				
                  %Pressure (kPa) = 10-1*(dBarr)
                  case 'PressureAbs', 
                     name = 'PRES';
                     data.PRES.values = (values{i})/10;  
					 data.PRES.comment = ['Pressure data in kPa converted to dBarr for toolbox'];	
			  
				  %Temperature (Celsius degree)
                  case 'Temperature', 
                    name = 'TEMP';
                    data.TEMp.values = values{1};
					data.TEMP.comment = ['Deg. Cel.'];
									  			  				  
				  %Conductivity (mS/cm) = 10-1*(S/m)
                  case 'Conductivity',
                      name = 'CNDC';
                      data.CNDC.values =(values{i})/10;
                      data.CNDC.comment = ['converted from mS/cm to S/m for Toolbox'];
					  
				  %DO (%)
                  %case 'AirSaturation(%)', 
                  %  name = 'DOXS';
                  %  data.DOXS.comment = ['Percent Dissolved Oxygen from O2 Concentration (uM)'];
                 
                 %% Doppler Current Meter - current Measurements similar to Nortek AquaProfiler - assumed no mag. decl. applied
                 
                 %Absolute Speed (cm/s) = 10-1*(m/s)
                  %case 'Abs Speed (cm/s)', 
                  %  name = 'CSPD';
                  %  data.CSPD.(values{i}) = data.(values{i})/10;
					%data.CSPD.comment = ['Absolute SeaWater Velocity converted from cm/s to m/s for Toolbox'];
					
				%Absolute Water Direction (Deg.M)
                 % case 'Direction(Deg.M)', 
                  %  name = 'CDIR';
					%data.CDIR.comment = ['Direction of the SeaWater Velocity'];
					
% 				%North(cm/s)  = 10-1*(m/s)
%                   case 'North(cm/s)', 
%                     name = 'VCUR';
%                     data.VCUR.(values{i}) = data.(values{i})/10;
% 					data.VCUR.comment = ['Northward water velocity converted from cm/s to m/ s for toolbox'];
% 					
% 				%East(cm/s)  = 10-1*(m/s)
%                   case 'East(cm/s)', 
%                     name = 'VCUR';
%                     data.UCUR.(values{i}) = data.(values{i})/10;
% 					data.UCUR.comment = ['Eastward water velocity converted from cm/s to m/ s for toolbox'];
% 					
% 				%Heading(Deg.M)
%                   case 'Heading(Deg.M)', 
%                     name = 'HEADING';
% 					data.HEADING.comment = ['Heading relative to Magnetic North'];
					
				%Tilt X(Deg)
                  case 'Pitch', 
                    name = 'PITCH';
					data.PITCH.comment = ['Horizontal movement'];
					
                 %Tilt Y(Deg)
                  case 'Roll', 
                    name = 'ROLL';
					data.ROLL.comment = ['Vertical movement'];
                 
                 %Signal Strength (backscatter - dB)
                  case 'PulseAttenuation', 
                    name = 'ABSI';
					data.ABSI.comment = ['Backscatter dB - Received Signal Strength'];
                    
                    %%%% there are also going to be ABSCI1,2,3,4? - 4 beam
                    %%%% unit which uses the best 3 beam solution
              
      end % end of Switch
  end % end of For loop
end % end of readData function
