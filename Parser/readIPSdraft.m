function sample_data = readIPSdraft(filename, mode )
%readIPSdraft Parses non-QC'd draft data exported from an IPS5 unit from ASL. 
%
% The estimated ice keel depth / draft is exported from the Ips5Extract software and
% stored in a 'IpsExtract_Output' folder within a .txt file with the calculated draft data.
%
% Variables Header for Draft export file is (space separated) : 
% YEAR MONTH DAY HOUR MINUTE SECOND DRAFT DERROR
%
% Filename example : ASL_IPS5_51104_Draft.txt
%
% Filename format : Company_Model_SerialNum_type(ice keel draft).txt
%
% this function also reads multiple *.csv sensor files from an IPS5 unit (ice profiling ADCP) 
% extracted via Ips5Extract v6.3.4 software, within the same folder as the *_Draft.txt file.
%
% The raw dataset is converted into 3 types of CSV files; range,sensor,and wave separated by phases 
% (up to 12 phases, though 5 to 8 is common for ArcticNet moorings).
%
% Filename example : iBO_BRG-17_201709_201910_p01_sensor_ed00.csv
%
% Filename format is : Program_MooringID_StartYearMonth_EndYearMonth_phase_type(range,sensor,wave)_ed00(version 00)
%
% The below descriptions are separated as such due to the timing of different types not synchronized amongst each type.
% Thus, the 3 types of data are separated out into 3 file types with which the data can be extracted from.
%
% RANGE CHANNEL DESCRIPTIONS
% 	Time [s]                  	- The logical record time specified as an offset in seconds from the start date.
% 	NumTarget                 	- The number of targets detected within a single ping. The instrument stores up to five targets per ping.
% 	Range:N [m]               	- The distance from the instrument acoustic transducer face to the Nth target. This is based on the sound speed defined by the user during instrument configuration.
% 	MaxAmplitude:N [counts]   	- The maximum echo amplitude achieved within the Nth target.
% 	Persistence:N [ms]        	- The timespan for which the target echo remained above the start and stop amplitude thresholds defined by the user during instrument configuration.
% 	MaxAmpIndRange [m]        	- The range corresponding to the maximum echo amplitude achieved within the full ping echo profile.
% 	MaxAmplitude [counts]     	- The maximum amplitude achieved within the full ping echo profile.
% 	BurstFlag                 	- Equal to 1 when the logical record corresponds to a ping acquired in burst mode. Otherwise, equal to 0.
% 	PingNumber [count]			- The unique instrument ping number used to correlate records in different files
%
% SENSOR CHANNEL DESCRIPTIONS
% 	SecondsSincePhaseStart [s] 	- The logical record time specified as an offset in seconds from the start date.
% 	TiltX [deg]                	- Beam tilt from vertical on the x-axis.
% 	TiltY [deg]                	- Beam tilt from vertical on the y-axis
% 	ParosTemperature [C]       	- Temperature from the Paros pressure sensor.
% 	Pressure [dbar]            	- Total external pressure.
% 	Temperature [C]            	- Ambient instrument temperature (optional).
% 	Battery [counts]           	- Battery voltage.
% 	PingNumber [count]			- The unique instrument ping number used to find corresponding records in different files.
%
% WAVE DESCRIPTIONS
% 	Time [s]             		- The logical record time specified as an offset in seconds from the start date.
% 	IndexInBurst         		- The relative sample index within a burst.
% 	BurstNum            		- The index of the burst.
% 	Amplitude [counts]   		- The maximum amplitude of the ping.
% 	MaxAmpIndex          		- The ping echo amplitude sample index corresponding to the maximum amplitude achieved within the ping echo profile. 
% 	Gain                 		- The gain used when during ping echo reception.
% 	Range [m]            		- The distance from the instrument acoustic transducer face to the first target. This is based on the sound speed defined by the user during instrument configuration.
% 	TiltX [deg]          		- Beam tilt from vertical on the x-axis.
% 	TiltY [deg]          		- Beam tilt from vertical on the y-axis
% 	Battery [V]          		- Battery voltage.
% 	Temperature [C]      		- Ambient instrument temperature (optional).
% 	ParosPressure [dbar] 		- Total external pressure.
% 	ParosTemperature [C] 		- Temperature from the Paros pressure sensor.
% 	Pressure [dbar]      		- Total external pressure.
% 	PingNumber [count]   		- The unique instrument ping number used to find corresponding records in different files.
%
%
% Inputs:
%   filename  	- Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 		 Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>			
%
% Copyright (c) 2020, Amundsen Science & ArcticNet
% http://www.amundsen.ulaval.ca/
% http://www.arcticnet.ulaval.ca/
% All rights reserved.
%
% Copyright (c) 2020, Australian Ocean Data Network (AODN) and Integrated 
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

[filepath,name,~] = fileparts(filename);
%extracting filename parts 
% Filename example : ASL_IPS5_51104_Draft.txt
unitInfo = textscan(name, '%s %s %d %s', 'Delimiter', '_');
%unitMake = char(unitInfo{1}); % company
%unitModel = char(unitInfo{2});

% extracting ice drafts - first approximation
data = readData(filename);
% extracting sensor data and combining the ice draft and sensor data tables
data = readSensorData(filepath,data);     
 
% copy all of the information over to the sample data struct
sample_data = struct;

sample_data.toolbox_input_file              = filename;
sample_data.meta.instrument_make            = char(unitInfo{1}); 
sample_data.meta.instrument_model           = char(unitInfo{2});
sample_data.meta.instrument_serial_no       = unitInfo{3};
%sample_data.meta.median_depth               = median(data.PRES.values);
sample_data.meta.instrument_sample_interval = median(diff(data.TIME.values*24*3600)); % varies by phase thus, not accurate for IPS
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


function data = readData(filename)
%READDATA Reads the sample data from the file.

  data = struct;
    
  fid = fopen(filename, 'rt');
  params = textscan(fid,'%s',8,'Delimiter',' ');
  values = textscan(fid,repmat('%f',[1,8]),'Delimiter',' ');
  fclose(fid);
  
%  YEAR MONTH DAY HOUR MINUTE SECOND DRAFT DERROR
% date and time columns 1-6 to make one TIME variable

Y=values{1};
M=values{2};
D=values{3};
H=values{4};
MN=values{5};
S=values{6};
Draft = values{7}; % draft values
Draft_err = values{8}; % Draft Error values

% adding DateTime to values and params
params{1,1}{9,1} = 'DATETIME';
values{9} = datenum(Y,M,D,H,MN,S); % DATETIME
DTime = values{9};

% Sorting the dateTime in ascending order without duplicates
[~,Sort_order]=unique(DTime,'sorted'); 

DTime     = DTime(Sort_order);
Draft     = Draft(Sort_order);
Draft_err = Draft_err(Sort_order);

%identify impossible dates, to not include in the final data struct
[c,~] = sort(DTime > datenum(2000,01,01));

DTime     = DTime(c);
Draft     = Draft(c);
Draft_err = Draft_err(c);

%sorting data struct by unique sorted time values
data.TIME.values=DTime;
data.ICE_DRAFT.values=Draft;
data.ICE_DRAFT_ERR.values=Draft_err;
%
% selecting +5m to -40m ice draft depth selection ; typically the Beaufort
% Sea is between 15-30m drafts. 60-100m is maybe more applicable for the 
% Labrador Sea. Thus, change to 60-100m if IPS are ever deployed in the 
% Eastern Arctic. - ShawnM Feb 19, 2020
%
inRangeA = Draft >=-5 & Draft <= 40; 

data.ICE_DRAFT.values = data.ICE_DRAFT.values(inRangeA)*(-1); %values inverted for vertical plotting sense
data.ICE_DRAFT.comment = 'Ice Keel Depth (m)';

data.ICE_DRAFT_ERR.values = data.ICE_DRAFT_ERR.values(inRangeA)*(-1); %values inverted for vertical plotting sense
data.ICE_DRAFT_ERR.comment = 'Ice Keel Depth Error (m)'; 

data.TIME.values = data.TIME.values(inRangeA); 
data.TIME.comment = 'DateTime';

% house cleaning 
clear Y M D H MN S params values Sort_order inRangeA DTime Draft Draft_err c

end % end of readData function

function data = readSensorData(filepath, data)
%% readIPS reads a directory with multiple.csv files from an IPS5 unit (ice profiling ADCP) 
% extracted via Ips5Extract v6.3.4 software.
%
% The raw dataset is converted into 3 types of CSV files; range,sensor,and wave separated by phases 
% (up to 12 phases, though 5-8 is common for ArcticNet moorings).
%
% Filename example : iBO_BRG-17_201709_201910_p01_sensor_ed00.csv
%
% Filename format is : Program_MooringID_StartYearMonth_EndYearMonth_type(range,sensor,wave)_ed00(version 00)
%
% The below descriptions are separated as such due to the timing of different types not synchronized amongst each type.
% Thus, the 3 types of data are separated out into 3 file types with which the data can be extracted from.
%
% SENSOR CHANNEL DESCRIPTIONS
% 	SecondsSincePhaseStart [s] 	- The logical record time specified as an offset in seconds from the start date.
% 	TiltX [deg]                	- Beam tilt from vertical on the x-axis.
% 	TiltY [deg]                	- Beam tilt from vertical on the y-axis
% 	ParosTemperature [C]       	- Temperature from the Paros pressure sensor.
% 	Pressure [dbar]            	- Total external pressure.
% 	Temperature [C]            	- Ambient instrument temperature (optional).
% 	Battery [counts]           	- Battery voltage.
% 	PingNumber [count]			- The unique instrument ping number used to find corresponding records in different files.
%
%% read all CSV files from the ice draft filepath
%  
myfiles =   dir(filepath);
filenamesCSV = {myfiles(:).name}; 
filefolders={myfiles(:).folder};
%  
%list only csv files with sensor in name
csvfiles=filenamesCSV(endsWith(filenamesCSV,'sensor_ed00.csv'));
csvfolders=filefolders(endsWith(filenamesCSV,'sensor_ed00.csv')); 
%
%Make a cell array of strings containing the full file locations of the
%files.
files=fullfile(csvfolders(:),csvfiles(:));
%
%Import the csv tables into out table
for i = 1:length(files)
    opts = detectImportOptions(files{i},'PreserveVariableNames', true,'VariableNamesLine',1);
    out{i} = readtable(files{i},opts);
    % maybe import col 1 applying a datenum conversion from str to ISO date
end
%
% merge the outputs into a usable table
values = vertcat(out{:});
%
% Getting variable names from the table
params = opts.VariableNames;  
nParams = length(params); 
%    
% Sorting the dateTime of the table in ascending order without duplicates
% new date table is created, though dates are strings 
%
isodates(:,1) = datenum(values{:,1},'yyyy/mm/dd-HH:MM:SS');
%
% new array values with mean of sorted dateTime
for j = 2: nParams
    isodates(:,j) = values{:,j};
end
%
% DateTimes from ICE_DRAFT and SENSOR data
T1 = datetime(data.TIME.values, 'convertFrom', 'datenum');
T2 = datetime(isodates(:,1), 'convertFrom', 'datenum');
%
% Create timetables of ICE_Draft and SENSOR data
TT2 = timetable(T2,isodates(:,2),isodates(:,3),isodates(:,4),isodates(:,5),isodates(:,7));
TT1 = timetable(T1, data.ICE_DRAFT.values, data.ICE_DRAFT_ERR.values );
%
% Union join based-on times of the 2 timetables binned by minute using mean
% while applying a linear interpolation between times.
%TT3 = synchronize(TT1,TT2,'union','mean');
TT3 = synchronize(TT1,TT2,'union');
TT4 = retime(TT3,'minutely','linear');
% 
% re-creating the final table for passing up for final IMOS/CF formating
t = timetable2table(TT4);
tt(:,1) = datenum(t{:,1});
for k = 2: width(t)
    tt(:,k) = t{:,k};
end
%
% finalizing the data struct for passing back up
data.TIME.values = tt(:,1);
data.ICE_DRAFT.values = tt(:,2);
data.ICE_DRAFT_ERR.values = tt(:,3);

data.ROLL.values = tt(:,4);
data.ROLL.comment = 'TiltX[deg]';

data.PITCH.values = tt(:,5);
data.PITCH.comment = 'TiltY[deg]';

data.TEMP.values = tt(:,6);
data.TEMP.comment = 'degrees celcius - at Paros ABS pressure sensor';

data.PRES.values = tt(:,7);
data.PRES.comment = 'Absolute external pressure (dbar)';

data.VOLT.values = tt(:,8);
data.VOLT.comment = 'Battery Voltage';  
%
%% House cleaning
clear values params nParams myfiles filenamesCSV filefolders csvfiles csvfolders files opts out i j k isodatesT1 T2 TT1 TT2 TT3 TT4 t tt

end % end of readData function