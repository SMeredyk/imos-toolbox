function sample_data = readMkvL( filename, mode )
% readMkvL reads a .RAW (CSV) data file retrieved from a JFE-ALEC Compact ACLW logger.
% via MDS software (ver 1.0.2 - 2002) - Only one channel of Light data.
%
% Inputs:
%   filename    - Cell string array passed from AlecParse.m, containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 		Shawn.Meredyk@arcticNet.ulaval.ca (ArcticNet - ULaval - Canada)
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
    rawText = textscan(fid, '%s', 'Delimiter', '');  % read all
    fclose(fid);  
    
    
    [header, channel, nChannels] = readHeader(rawText{1}); % read Header
    data = readData(filename, channel, nChannels);  % read Data and convert
    
catch e
    if fid ~= -1, fclose(fid); end
    rethrow(e);
end

% copy all of the information over to the sample data struct
sample_data = struct;

sample_data.toolbox_input_file              = filename;
sample_data.meta.instrument_make            = 'JFE_ALEC';
sample_data.meta.instrument_model           = header.InstType;
sample_data.meta.instrument_serial_no       = header.InstNo;
sample_data.meta.instrument_sample_interval = median(diff(data.TIME.values*24*3600));
%sample_data.meta.instrument_burst_interval  = header.BurstTime; % seconds between bursts
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

function [header, channel,nChannels] = readHeader(rawText)
%READHEADER Reads the header and coefficents from the file. 
  header = struct;
      
  startHeader   = '[Head]';
  endHeader     = '[Coef]'; %all file types have [Item] except compactActw , which has [Data]
  fmtHeader     = '%s%s';
  delimHeader   = '=';
  
  iStartHeader = find(strcmp(startHeader, rawText)) + 1;
  iEndHeader = find(strcmp(endHeader, rawText))-1;
    
  headerCell = rawText(iStartHeader:iEndHeader);
  nFields = length(headerCell);
  for i=1:nFields
      tuple = textscan(headerCell{i}, fmtHeader, 'Delimiter', delimHeader);
      if 1 == isempty(tuple{2}), tuple{2} = {'NAN'}; end
        
      header.(tuple{1}{1}) = tuple{2}{1};
  end   % end of for loop

%Now the header function reads the coefficients, to create the channel struct

  channel	= struct;
    
  startCoefficient 	= '[Coef]';
  endCoefficient 	= '[Item]';
  delimCoef         = ',';
  fmtCoef           = '%f%f%f%f'; % 4 coefficients
  
  iStartCoef 	= find(strcmp(startCoefficient, rawText)) + 1;
  iEndCoef      = find(strcmp(endCoefficient, rawText))-1;
 
  
  coefCell 	= rawText(iStartCoef:iEndCoef);
  nChannels = length(coefCell);     % number of channels
  
  for i=1:nChannels
      coef 	= textscan(coefCell{i}, fmtCoef,'Delimiter', delimCoef);
      
      for j=1:4
          eval(['channel.Ch' num2str(i) '(j) = coef{j};']); end 
  end   %end for loop
  
end % end of readHeader function


function data = readData(filename, channel, nChannels)
%READDATA Reads the sample data from the file.

  data = struct;
  
  dataDelim = ',';
  dataFmt = ['%s', repmat('%f', 1, nChannels)];
   
  fid = fopen( filename, 'rt' );
    c = 1;
    headerL{c} = fgetl(fid); 
    
    % Read file until [Item]   
    while isempty( strfind( headerL{c}, '[Item]'))
        c = c + 1; 
        headerL{c}= fgetl(fid);
    end
    
    % Read the data
    values = textscan( fid, dataFmt, 'delimiter', dataDelim );
    fclose( fid );

    % convert the raw data into real values      
    data.TIME.values = datenum(values{1},'yyyy/mm/dd HH:MM:SS');
    data.TIME.comment = ['Time format - yyyy/mm/dd HH:MM:SS'];
    
    data.PAR.values = channel.Ch1(2).*(values{2});
    data.PAR.comment = ['Light Quantum : [umol_phtn/m2/s]']; 
    		  
end % end of readData function
