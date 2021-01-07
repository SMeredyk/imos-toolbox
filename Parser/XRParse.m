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
% Author : 		 Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributors : Guillaume Galibert <guillaume.galibert@utas.edu.au>			
%
% Copyright (c) 2017, Amundsen Science & ArcticNet
% http://www.amundsen.ulaval.ca/
% http://www.arcticnet.ulaval.ca/
% All rights reserved.
%
% Copyright (C) 2017, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.
% If not, see <https://www.gnu.org/licenses/gpl-3.0.en.html>.
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

fname = filename; % to accomodate RSKTools functions without having to varName replace - shawn Jan 6, 2021

if strcmpi(ext, '.dat') && strcmp(line(1:9), 'Model=RBR')
    % use the classic XR420 parser for RBR Windows v 6.13 file format
    sample_data = readXR420(filename, mode);
elseif strcmpi(ext, '.txt') && strcmp(line(1:12), 'Model=XR-420')
            sample_data = readXR620(filename, mode); % This is a newer Ruskin exported XR420 dataset, 
			% imports CT, CT-FL-Tu-DO ; no salinity or depth imported
elseif strcmpi(ext, '.txt') && strcmp(line(1:12), 'Model=XR-620')
            sample_data = readXR620(filename, mode); % This is a newer Ruskin exported XR420 dataset, 
			% imports CT, CT-FL-Tu-DO ; no salinity or depth imported		
elseif strcmpi(ext, '.txt') && strcmp(line(1:12), 'Model=RBRcon')		
		   	sample_data = readXRConcertoDuo(filename, mode); % This is a Concerto model 
			%this parser is the same as XR620
elseif strcmpi(ext, '.txt') && strcmp(line(1:10), 'Model=RBRd')			
		   	sample_data = readXRConcertoDuo(filename, mode); % This is a Duo model - 
			%this parser is the same as XR620
else strcmpi(ext, '.rsk')			
		   	sample_data = readRSKfile(fname, mode); % This is any RBR unit with data in .rsk format
			%this is a very different / simple parser in comparison to the other model type functions.
			% this may become the default parser as legacy exporting is limited to legacy products.
end % end of main
