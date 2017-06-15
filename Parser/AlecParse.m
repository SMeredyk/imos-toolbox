function sample_data = AlecParse(filename, mode)
% AlecPARSE reads a data file retrieved from a variety of JFE-ALEC instrument
% loggers.
%
% This function is able to read in a single file retrieved from an JFE-ALEC
% ACTW-CMP, ACTW-USB, ALW-CMP, ACLW-CMP,ACLW-USB, ACT-HR generated using 
% ALEC software and direct the toolbox to the respective parser. 
% The pressure data is returned in a sample_data struct.
%
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 	    Shawn.Meredyk@arcticNet.ulaval.ca (ArcticNet - ULaval - Canada)
% Contributors:	Guillaume Galibert <guillaume.galibert@utas.edu.au>, 
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

% ensure that there is exactly one argument, 
% and that it is a cell array of strings
narginchk(1,2);

if ~iscellstr(filename), error('filename must be a cell array of strings'); end

% only one file supported currently
filename = filename{1};
if ~ischar(filename), error('filename must contain a string'); end

% open the file, and read in the data
try
    fid = fopen(filename, 'rt');
    rawText = textscan(fid, '%s', 'Delimiter', '');
    fclose(fid);
    
    header = readHeader(rawText{1});
    
catch e
    if fid ~= -1, fclose(fid); end
    rethrow(e);
end

% copy all of the information over to the sample data struct
sample_data = struct;

sample_data.toolbox_input_file              = filename;
sample_data.meta.instrument_make            = 'JFE_ALEC';
sample_data.meta.featureType                = mode;

%%%% Trying to assign model and serial_no from header data from a variety of file possibilities
%%%% There are three main file types that either use SondeName or SensorType or InstType
%%%% The full header needs to be read as all file types have [Item] except compactActw , which has [Data]

% to set filename back to a cell string array, for passing to next parser.

if isfield(header,'SondeName')== 1  % If the header scan filled-in the variable SondeName, otherwise elseif
% Infinity Series Units
sample_data.meta.instrument_model           = header.SondeName;
sample_data.meta.instrument_serial_no       = header.SondeNo;
  
	% Infinity Model and Serial number
	equipname = regexp(header.SensorType,'T0K0U3B0');
    if equipname == 1
	    sample_data = readInfinityAclw(filename, mode); % This is an ACLW-USB model 
	else
	    sample_data = readInfinityActw(filename, mode); % This is an ACTW-USB model 'T1C1C2B0' 
    end  
  
elseif isfield(header, 'SensorType')== 1 % if the header scan filled-in variable SensorType, otherwise elseif
% Compact Series Units (Compact-CTW Software V1.02 - 2009?)
sample_data.meta.instrument_model           = header.SensorType;
sample_data.meta.instrument_serial_no       = header.SerialNo;

	% Compact Model and Serial number
	equipname = regexp(header.SensorType,'TCWWWE');
	if equipname ==1
		sample_data = readCompactActw(filename, mode); % This is a Compact-CTW model
	else 
		header.SensorType = 'TKURB';
		sample_data = readCompactAclw(filename, mode); % This is a Compact-CLW model
	end


elseif isfield(header,'InstType') == 1
            sample_data.meta.instrument_model       = header.InstType;
            sample_data.meta.instrument_serial_no   = header.InstNo;
            
	switch header.InstType
               
		case 'QB'
            sample_data = readCompactAlw(filename, mode); % This is an ALW-USB model
		
		case 'TC'
			sample_data = readCompactActHr(filename, mode); % This is an ACT-HR model 'TC' 
			
		case 'TKURB'
           	sample_data = readCompactAclw(filename, mode); % This is a Compact-CLW model

	otherwise, error('unknown type'); 
end 	% end of switch

end 	% end of InstType
	


function header = readHeader(rawText)
%READHEADER Reads the header from the file.
  header = struct;
      
  startHeader   = '[Head]';
  endHeader     = {'[Coef]' '[Item]' '[Data]'}; %all file types have [Item] except compactActw , which has [Data]
  fmtHeader     = '%s%s';
  delimHeader   = '=';
  
  %iStartHeader = find(strcmp(startHeader, rawText)) + 1;
  iStartHeader = find(cellfun('isempty', regexp(rawText, '^\[Head]|') ) == 0 ) + 1;
  iEndHeader = find(cellfun('isempty', regexp(rawText, '^\[Coef]|^\[Item\]|^\[Data\]') ) == 0 ) - 1;   
      
  headerCell = rawText(iStartHeader:iEndHeader(1));
  nFields = length(headerCell);
  for i=1:nFields
      tuple = textscan(headerCell{i}, fmtHeader, 'Delimiter', delimHeader);
      
      if 1 == isempty(tuple{2}), tuple{2} = {'NAN'}; end
      
      tuple{1}{1}=strrep(tuple{1}{1},' ','_');  %% Added Aforest 30-Jan-2107 for valid field name (no space)
      header.(tuple{1}{1}) = tuple{2}{1};
     
  end
end     % readHeader function end

end % closes main function