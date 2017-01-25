function sample_data = XRParse( filename, mode )
%XRPARSE Parses a RBR data file retrieved from an RBR XR420 (CT, CT-FL-Tu-DO), 
% Concerto, Duo or XR620 data logger.
%
% This function is able to read in a single file retrieved from an RBR
% XR420, Concerto/Duo or XR620 data logger generated using RBR Windows v 6.13 software 
% or Ruskin software. The pressure data is returned in a sample_data
% struct. The salinty data and Speed of Sound are not imported for the XR420 or 
% XR620 models; only Concerto/Duo because of the onboard CTD data.
%
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode.
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 		Guillaume Galibert <guillaume.galibert@utas.edu.au>
% Contributors :Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
%
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

% ensure that there is exactly one argument, 
% and that it is a cell array of strings
narginchk(1,2);

if ~iscellstr(filename), error('filename must be a cell array of strings'); end

% only one file supported currently
filename = filename{1};
if ~ischar(filename), error('filename must contain a string'); end

[~, ~, ext] = fileparts(filename);

% read first line in the file
line = '';
try
    fid = fopen(filename, 'rt');
    line = fgetl(fid);
    fclose(fid);
    
catch e
    if fid ~= -1, fclose(fid); end
    rethrow(e);
end


if strcmpi(ext, '.dat') && strcmp(line(1:9), 'Model=RBR')
    % use the classic XR420 parser for RBR Windows v 6.13 file format
    sample_data = readXR420(filename, mode);
elseif strcmpi(ext, '.txt') && strcmp(line(1:12), 'Model=XR-420')
            sample_data = readXR620(filename, mode); % This is a neewer Ruskin exported XR420 dataset, 
			% imports CT, CT-FL-Tu-DO ; no salinity or depth imported
elseif strcmpi(ext, '.txt') && strcmp(line(1:12), 'Model=XR-620')
            sample_data = readXR620(filename, mode); % This is a neewer Ruskin exported XR420 dataset, 
			% imports CT, CT-FL-Tu-DO ; no salinity or depth imported		
elseif strcmpi(ext, '.txt') && strcmp(line(1:12), 'Model=RBRcon')		
		   	sample_data = readXRConcertoDuo(filename, mode); % This is a Concerto model 
			%this parser is the same as XR620
else strcmpi(ext, '.txt') && strcmp(line(1:10), 'Model=RBRd')			
		   	sample_data = readXRConcertoDuo(filename, mode); % This is a Duo model - 
			%this parser is the same as XR620

end % end of main
