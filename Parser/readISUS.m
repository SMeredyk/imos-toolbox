function sample_data = readISUS( foldername, mode )
%readISUS Parses several data files retrieved from a Satlantic nitrate sensor ISUS V2 /V3 
%
% This function is able to read a folder of .dat files retrieved from an ISUS unit.
% The data is recovered via the ISUSCom 2.1.6 software
%
% Inputs:
%   foldername  - Cell array containing the name of the folder to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
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

if ~ischar(foldername), error('foldername must contain a string'); end
 

 
%read folder contents : directory, filepattern, filterout, fpattern
%directory = foldername;
filepattern = 'DAT';
filterout = 0;
fpattern = 'DIVE';

[n, datafiles] = get_folder_content(foldername, filepattern, filterout, fpattern);

%now n is number of file names listed in structure files
% to access file names use this:

%datafiles(1).name % this will give you just the file name of the first file in the list
%datafiles(1).fullname %this will give you absolute file name ... so with absolute path
 
%read the various DAT files to create a consolidated data structure 
data = readData(n, datafiles);  % read Data and convert

% copy all of the information over to the sample data struct
sample_data = struct;

sample_data.toolbox_input_file              = datafiles(1).fullname;
sample_data.meta.instrument_make            = 'Satlantic';
sample_data.meta.instrument_model           = 'V3';
sample_data.meta.instrument_serial_no       = '135'; % need to extract serial number from main folder name
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

function data = readData(n, datafiles)
%READDATA Reads the sample data from the numerous files in the specified folder.

  data = struct;
  
  dataDelim = ',';
  dataFmt = '%*s %s %f %f %*[^\n]'; % from pascal`s code
  metahdr = 'SATNHR';
  
  for p=1:n
  
    % Read the data
	fid = fopen( datafiles(p).name, 'rt' );
    values = textscan( fid, dataFmt, 'delimiter', dataDelim, 'commentStyle',metahdr ); % modified from pascal`s code
    fclose( fid );

	% Building Date and Time from DAT file
	date = values{ 1 };
	date = char( date );
	doy = str2num( date(:, 5: end) );
	date = str2num( date(:, 1: 4) );
	date = datestr( doy2date( doy, date ) );

	%hours and minutes
	hours = values{ 2 };
	hh = floor( hours );
	minu = rem(hours,hh) * 60;

	% ss = ones( size( minu ) );
	I =( floor( minu ) == 0 );
	ss( I ) = minu( I ) * 60;
	ss( ~I ) = rem( minu(~I), floor( minu(~I) ) ) * 60;
	tmp = datevec( date );
	tmp(:,4) = hh;
	tmp(:,5) = floor(minu);
	tmp(:,6) = ss;
	tmp( isnan( tmp ) ) = 0;
	date = datestr( datenum( tmp ), 'dd-mmm-yyyy HH:MM:SS.FFF' );

	% importing nitrate data into values struct
	nitrates = values{3};
	
	  % Remove values == 0;
		date( nitrates == 0, : ) = [];
		nitrates( nitrates == 0 ) = [];
		
	%re-import date and nitrates to value struct after zeros removed
	values{1} = date;
	values{3} = nitrates;
	
    % convert the raw data into real values      
    data.TIME.values = datenum(values{1});
    data.TIME.comment = 'Time'; 
     
    data.NTRA.values = values{3};
    data.NTRA.comment = ['mole_concentration_of_nitrate_in_sea_water : [mole l-1]']; 
	
  end            

   			  
end % end of readData function


% Function return the folder content filtered
% function [ n, datafiles ] = get_folder_content( directory, filepattern, filterout, fpattern)
% 
% input:
%        - directory
%        - filepattern 
%        - filterout    |  1 or 0
%        - fpattern     |  not include files with this pattern
%          fpattern is charracter array created by fpattern =
%          char('string1', 'string2', ...,'stringn');
%          or in case of one string fpattern = 'string';
%        
% output:
%        - n | number of files
%        - datafiles | structure with datafiles
%          datafiles(i).name -- sigle file name
%          datafiles(i).fullname -- including full path 

function [ n, datafiles ] = get_folder_content( foldername, filepattern, filterout, fpattern)

%initialisation of the structure
datafiles(1).name = '';
datafiles(1).fullname = '';
  
  if nargin <3,
     filterout = 0; 
  end
  
     tmpdir =  dir([foldername filepattern '*']);
     nn = length(tmpdir);
     %filtering not necessary files
     cnt = 0;
  display('Folder contain files for further processing:')   
  if filterout==1
     [nstr s] = size(fpattern);
     filter_flag = 0;
     for i=1:nn
         for st=1:nstr
            S = regexp(tmpdir(i).name,fpattern(st,:));
            if(isempty(S))
                filter_flag=0;
            else
                filter_flag=1;
            end
         end
         if(filter_flag==0)
             % filter pattern not maching ... so this is correct file
             cnt = cnt+1;
             datafiles(cnt).name = tmpdir(i).name;
             datafiles(cnt).fullname = [foldername tmpdir(i).name];
             
         else
             %filter pattern matching 
             %display(['Filtering out: ' tmpdir(i).name ' file'])
         end
     end
     
     
  else
     for i=1:nn
         datafiles(i).name = tmpdir(i).name;
         datafiles(i).fullname = [foldername tmpdir(i).name];
     end
      
  end

  
  n = length(datafiles);
  for i=1:n
      display([ num2str(i) ' file: ' datafiles(i).name ])
  end
  display(['Number of listed files: ' num2str(n)]);
  return

end  


