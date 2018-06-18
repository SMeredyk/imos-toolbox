function sample_data = AanderaaParse( filename, mode )
% AanderaaPARSE reads a tab delimited data file converted from the 5059 Data 
% reading program , for RCM4,7,11 and seaguard data loggers.
%
% This function is able to read in a single file retrieved from Aanderaa
% RCM4,7,11 and Seaguard data loggers generated using 
% 5059 Data Reading Program and directs the toolbox to either the 'rcm' or 
% 'seaguard' parsers. The pressure data is returned in a sample_data struct.
%
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author :      Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributors: Guillaume Galibert <guillaume.galibert@utas.edu.au>
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

if strcmpi(ext, '.txt') && strcmp(line(1:7), 'Model:S')
    % use the Seaguard file format
    sample_data = readSeaguard(filename, mode);
elseif strcmpi(ext, '.txt') && strcmp(line(1:9), 'Model:RCM')
    % use the RCM4,7,11 file format
    sample_data = readRCM(filename, mode);
else 
    % use the RDCP file format
    sample_data = readRDCP(filename, mode);
end

