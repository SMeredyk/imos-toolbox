function sample_data = readRDCP(filename, mode)
% readRDCP reads a data file (*.txt) retrieved from a Seaguard RDCP data
% logger.RDCP software flags such as {3} are removed along with the record
% number and interval columns (exported by RDCP Studio software).
%
% The pressure data is returned in a sample_data struct. 
% 
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Sample Header Format (below):
% Model::RDCP
% SerialNumber::270
% RefReading::593
% BlankDist::1
% CellSize::4
% NumCells::30
% CellOverlap::50
% DistToFirstBin::3
% RecordInterval::3600
% StartTime::14.08.2007_01:00
% 
% [Data]
%
% Author :      Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributors: Guillaume Galibert <guillaume.galibert@utas.edu.au>
%               Pascal_Guillot@uqar.ca (UQAR - Canada)			
%
% Copyright (c) 2017, eMarine Information Infrastructure (eMII) and Integrated 
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

% Correction for pressure offset in air - Originally added by AForest 27-Jan-2017 with
% comments for history on 30-Jan-2017 to Nortek current profiler toolbox code.

% based on first 5 measurements within 10 m range
[~,NAME,~] = fileparts(filename);
first_mes=data.PRES.values(1:5);
first_mes=first_mes(first_mes<10);

if  ~isnan(first_mes)
    disp(['Please note: ', NAME,': pressure offset in air : ',...
        num2str(ceil(max(first_mes))),'-dbar Pressure Offset Applied']);
    
	%pressure=pressure-mean(first_mes);
    data.PRES.values=data.PRES.values-mean(first_mes);
    
	% Commenting the Metadata history
    PressureOffsetComment=[mfilename,'.m: Raw pressure data from ', NAME,...
        ' was corrected for a pressure offset in air of ',...
        num2str(round(mean(first_mes),1)),'dbar'];
    
    sample_data.history = sprintf('%s - %s', ...
            datestr(now_utc, readProperty('exportNetCDF.dateFormat')), ...
            PressureOffsetComment);
else
    disp(['Please note: ', NAME,': pressure offset in air : ',...
        num2str(ceil(max(data.PRES.values(1:5)))),...
        '-dbar and NO pressure offset was applied']);
end

sample_data.toolbox_input_file              = filename;
sample_data.meta.instrument_make            = 'Aanderaa';
sample_data.meta.instrument_model           = 'RDCP';
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
    %sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc([data.(fields{k}).values{:}]);
    sample_data.variables{end}.coordinates  = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
    sample_data.variables{end}.comment      = data.(fields{k}).comment;
end	
end %end of main function

function [header, iData] = readHeader(rawText)
%READHEADER Reads the header , including channel coefficients from the file.

  header 	= struct;
  iData 	= [];
  
  startHeader 	= 'Model::RDCP';	
  endHeader 	= '[Data]';    
  fmtHeader  	= '%s%s';
  delimHeader 	= '::';
  
  
  iStartHeader = find(strcmp(startHeader, rawText));
  iEndHeader = find(strcmp(endHeader, rawText))-1;
  iData = iEndHeader + 2;  %data headers - 2 line after [Data]
  
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
  nParams = length(iParams{1})+1; % needs to see one other field?
  paramsFmt = repmat('%s', 1, nParams);
  params = textscan(params{1}, paramsFmt, 'Delimiter', dataDelim);
  dataFmt = ['%s', repmat('%f', 1, nParams-1)];
  values = textscan(fid, dataFmt, 'Delimiter', dataDelim);
  fclose(fid);
  
  
  for i=1:nParams
      switch params{i}{1} 
                                  
                 %Date Time (dd.mm.yy  HH:MM)
                  case 'Time Tag (Gmt)'
                    data.TIME.values = datenum(values{i}, 'dd.mm.yyyy  HH:MM');
                    data.TIME.comment = ['dd.mm.yyyy  HH:MM'];		
                    
                 %Reference Parameter (unitless) - not used by Toolbox
                  %case 'Reference' 
                  %  name = 'REF';
                  %  data.REF.values = values{i};
                  %  data.REF.comment = ['calibration reference value from Aanderaa'];
                    
                  %Battery (Volts)
                  case 'BatteryVoltage' 
                    name = 'VOLT';
                    data.VOLT.values = values{i};
					data.VOLT.comment = ['Volts'];
						
                  %Pressure (kPa) = 10-1*(dBarr)
                  case 'PressureAbs' 
                     name = 'PRES';
                     data.PRES.values = (values{i})/10;  
				     data.PRES.comment = ['Pressure data in kPa converted to dBarr for toolbox'];	
                     
                  %Pressure (kPa) = 10-1*(dBarr)
                  case 'PressureRel' 
                     name = 'PRES_REL';
                     data.PRES_REL.values = (values{i})/10;  
			         data.PRES_REL.comment = ['Pressure data in kPa converted to dBarr for toolbox'];
                     
				  %Temperature (Celsius degree)
                  case 'Temperature' 
                    name = 'TEMP';
                    data.TEMP.values = values{i};
					data.TEMP.comment = ['Deg. Cel.'];
									  			  				  
				  %Conductivity (mS/cm) = 10-1*(S/m)
                  case 'Conductivity'
                      name = 'CNDC';
                      data.CNDC.values =(values{i})/10;
                      data.CNDC.comment = ['converted from mS/cm to S/m for Toolbox'];
					                  
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
 				%Heading (deg.)
                   case 'Heading' 
                     name = 'HEADING';
                     data.HEADING.values = values{i};
 					 data.HEADING.comment = ['Magnetic Compass Heading'];
					
				%Tilt X(Deg)
                  case 'Pitch' 
                    name = 'PITCH';
                    data.PITCH.values = values{i};
					data.PITCH.comment = ['Horizontal movement'];
					
                 %Tilt Y(Deg)
                  case 'Roll' 
                    name = 'ROLL';
                    data.ROLL.values = values{i};
					data.ROLL.comment = ['Vertical movement'];
                 
                 %Signal Strength (backscatter - dB)
                  case 'PulseAttenuation' 
                    name = 'ABSI';
                    data.ABSI.values = values{i};
					data.ABSI.comment = ['Backscatter dB - Received Signal Strength'];
                    
                    %%%% ShawnM - Jun 9 - 2017
                    %%%% The dataset contains the speed and direction of
                    %%%% each beam at each bin depth. Code to extract this
                    %%%% information , will go here eventually.
              
      end % end of Switch
  end % end of For loop
end % end of readData function
