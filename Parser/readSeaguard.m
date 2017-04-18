function sample_data = readSeaguard(filename, mode)
% readSeaguard reads a data file retrieved from a Seaguard data logger.
%
% This function is able to read in a single text file (comma delimited) 
% retrieved from a Seaguard data logger exported (tab delimited) by the 
% Seaguard Studio software. Tab delimited filew needs to be converted into
% commas prior to running toolbox. The pressure data is returned in a sample_data struct. 
% 
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Header Example:
% Model:Seaguard
% ProductNumber:4430 
% SerialNumber:30
%
% Author :      Guillaume Galibert <guillaume.galibert@utas.edu.au>
% Contributors: Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
%               Pascal_Guillot@uqar.ca (UQAR - Canada)			
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

% Correction for pressure offset in air - Added AForest 27-Jan-2017 with
% comments for history on 30-Jan-2017

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
sample_data.meta.instrument_model           = header.Model;
sample_data.meta.instrument_serial_no       = header.SerialNumber;
sample_data.meta.head                       = 2000;     %1900 - 2000 kHz ; 25W in 1ms pulses ; 2 sec 
sample_data.meta.instrument_sample_interval = median(diff(data.TIME.values*24*3600));
sample_data.meta.beam_angle                 = 2;    % http://www.aanderaa.com/media/pdfs/zpulse-doppler-current-sensor-seaguard.pdf
sample_data.meta.binSize                    = 1.8; 
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
  
  startHeader 	= 'Model:Seaguard';	
  endHeader 	= '[Data]';    
  fmtHeader  	= '%s%s';
  delimHeader 	= ':';
  
  
  iStartHeader 	= find(strcmp(startHeader, rawText));
  iEndHeader 	= find(strcmp(endHeader, rawText))-1;
  iData 		= iEndHeader + 2;   %data headers-2 lines after [Data]
  
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
  nParams = length(iParams{1})+1; % needs to see one other field
  paramsFmt = repmat('%s', 1, nParams);
  params = textscan(params{1}, paramsFmt, 'Delimiter', dataDelim);
  dataFmt = ['%s', repmat('%f', 1, nParams-1)];
  values = textscan(fid, dataFmt, 'Delimiter', dataDelim);
  fclose(fid);
  
  for i=1:nParams
      switch params{i}{1}
                  %Date Time (dd.mm.yyyy HH:MM:SS)
                  case 'Time tag (Gmt)'
                    name = 'TIME';
                    data.TIME.values = datenum(values{i}, 'dd.mm.yyyy HH:MM:SS');
                    data.TIME.comment = ['dd.mm.yyyy HH:MM:SS'];		
              
                  %Battery (Volts)
                  case 'Battery Voltage (V)', 
                    name = 'VOLT';
                    data.VOLT.values = values{i};
					data.VOLT.comment = [''];
				
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
					data.TURB.comment = ['Turbidity data '...
					  'computed from bio-optical sensor raw counts measurements. The '...
					  'turbidity sensor is equipped with a 880nm peak wavelength LED to irradiate and a '...
					  'photodetector paired with an optical filter which measures everything '...
					  'that backscatters in the region of 650nm to 1000nm.']';
                  
				  % Two pressure resolution options MPa or kPa, depending on .cdb file
                  % export from 5059 software
                  
				  %Pressure (MPa) = 100-1*(dBarr)
                  case 'Pressure(MPa)', 
                     name = 'PRES';
                     data.PRES.values = (values{i})/100;   
					 data.PRES.comment = ['Pressure data converted from MPa to dBarr for toolbox'];	
				
                  %Pressure (kPa) = 10-1*(dBarr) ?
                  case 'Pressure(kPa)', 
                     name = 'PRES';
                     data.PRES.values = (values{i})/10;   
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
 
				  %DO (%)
                  case 'AirSaturation(%)', 
                    name = 'DOXS';
                    data.DOXS.values = values{i};
                    data.DOXS.comment = ['Percent Dissolved Oxygen from O2 Concentration (uM)'];
                 
                 %% Single-Point Current Meter - Basic current Measurements - assumed no mag. decl. applied
                 
                 %Absolute Speed (cm/s) = 100-1*(m/s)
                  case 'Abs Speed(cm/s)', 
                    name = 'CSPD';
                    data.CSPD.values = (values{i})/100; 
					data.CSPD.comment = ['Absolute SeaWater Velocity m/s'];
					
				%Absolute Water Direction (Deg.M)
                  case 'Direction(Deg.M)', 
                    name = 'CDIR';
                    data.CDIR.values = values{i};
					data.CDIR.comment = ['Direction of the SeaWater Velocity'];
					
				%North(cm/s)  = 100-1*(m/s)
                  case 'North(cm/s)', 
                    name = 'VCUR';
                    data.VCUR.values = (values{i})/100;
					data.VCUR.comment = ['Northward water velocity converted from cm/s to m/ s for toolbox'];
					
				%East(cm/s)  = 10-1*(m/s)
                  case 'East(cm/s)', 
                    name = 'UCUR';
                    data.UCUR.values = (values{i})/100;
					data.UCUR.comment = ['Eastward water velocity converted from cm/s to m/ s for toolbox'];
					
				%Heading(Deg.M)
                  case 'Heading(Deg.M)', 
                    name = 'HEADING';
                    data.HEADING.values = values{i};
					data.HEADING.comment = ['Heading relative to Magnetic North'];
					
				%Tilt X(Deg)
                  case 'Tilt X(Deg)', 
                    name = 'PITCH';
                    data.PITCH.values = values{i};
					data.PITCH.comment = ['Horizontal movement'];
					
                 %Tilt Y(Deg)
                  case 'Tilt Y(Deg)', 
                    name = 'ROLL';
                    data.ROLL.values = values{i};
					data.ROLL.comment = ['Vertical movement'];
                 
                 %Signal Strength (backscatter - dB)
                  case 'Strength(dB)', 
                    name = 'ABSI';
                    data.ABSI.values = values{i};
					data.ABSI.comment = ['Backscatter dB - Received Signal Strength'];
              
      end % end of Switch
  end % end of For loop
  


end % end of readData function
