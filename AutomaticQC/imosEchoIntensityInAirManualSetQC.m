function [sample_data, varChecked, paramsLog] = imosEchoIntensityInAirManualSetQC( sample_data, auto )
% IMOSECHOINTENISTYINAIRMANUALSETQC Quality control procedure for Teledyne Workhorse (and similar)
%
% Quality control ADCP instrument data, flaging velocity data based on echo
% intensity in air (if available)
%
% The user is prompted to select a given echo intensity threshold value based on
% deployment statistics so data with low signal to noise can be flagged as bad.
% Some knowledge on underwater accoustic of a given environementis advisable.
%
% Inputs:
%   sample_data - 	struct containing the entire data set and dimension data.
%   auto 		- 	logical, run QC in batch mode
%
% Outputs:
%   sample_data - 	same as input, with QC flags added for variable/dimension
%                 	data.
%   varChecked  - 	cell array of variables' name which have been checked
%   paramsLog   - 	string containing details about params' procedure to include in QC log
%
% Author:       Alexandre Forest <Alexandre.Forest@arcticnet.ulaval.ca>
% Contributor:  Guillaume Galibert <guillaume.galibert@utas.edu.au>,
%				Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
%
% Copyright (c) 2019, Amundsen Science & ArcticNet
% http://www.amundsen.ulaval.ca/
% http://www.arcticnet.ulaval.ca/
% All rights reserved.
%
% Copyright (c) 2019, Australian Ocean Data Network (AODN) and Integrated
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
narginchk(1, 2);
if ~isstruct(sample_data), error('sample_data must be a struct'); end

% auto logical in input to enable running under batch processing
if nargin<2, auto=false; end

varChecked = {};
paramsLog  = [];

% get all necessary dimensions and variables id in sample_data struct
idUcur = 0;
idVcur = 0;
idWcur = 0;
idCspd = 0;
idCdir = 0;
idABSIC = cell(4, 1);
for j=1:4
    idABSIC{j}  = 0;
end
lenVar = length(sample_data.variables);
for i=1:lenVar
    paramName = sample_data.variables{i}.name;
    
    if strncmpi(paramName, 'UCUR', 4),  idUcur = i; end
    if strncmpi(paramName, 'VCUR', 4),  idVcur = i; end
    if strcmpi(paramName, 'WCUR'),      idWcur = i; end
    if strcmpi(paramName, 'CSPD'),      idCspd = i; end
    if strncmpi(paramName, 'CDIR', 4),  idCdir = i; end
    for j=1:4
        cc = int2str(j);
        if strcmpi(paramName, ['ABSIC' cc]), idABSIC{j} = i; end
    end
end

% Also get time % AForest 16 February 2017
tTime = 'dimensions';
iTime = getVar(sample_data.(tTime), 'TIME');
if iTime == 0
    tTime = 'variables';
    iTime = getVar(sample_data.(tTime), 'TIME');
    if iTime == 0, return; end
end
time = sample_data.(tTime){iTime}.data;

% check if the data is compatible with the QC algorithm
idMandatory = iTime & (idUcur | idVcur | idWcur | idCspd | idCdir);

% check if we have nortek
if idABSIC{j}==0
    beams=3; %Nortek
else
    beams=4; %RDI
end

for j=1:beams
    idMandatory = idMandatory & idABSIC{j};
end
if ~idMandatory, return; end

qcSet           = str2double(readProperty('toolbox.qc_set'));
badFlag         = imosQCFlag('bad',             qcSet, 'flag');
goodFlag        = imosQCFlag('good',            qcSet, 'flag');
rawFlag         = imosQCFlag('raw',             qcSet, 'flag');

%Pull out echo intensity
sizeData = size(sample_data.variables{idABSIC{1}}.data);
ea = nan(beams, sizeData(1), sizeData(2));
for j=1:beams
    ea(j, :, :) = sample_data.variables{idABSIC{j}}.data;
end

%% Get the pre-existing NoiseFloor data from the DDB
noiseFloor      = 20; % assuming standard 20 counts , if no data in DB
%nSample         = length(sample_data);

% get badbins for each data set
%for k = 1:nSample
      
    %check to see if noiseFloor data are available already from the ddb
    %if isfield(sample_data{k}.meta, 'deployment')
        if ~isempty(sample_data.meta.deployment.NoiseFloor)
            noiseFloor = sample_data.meta.deployment.NoiseFloor;
        end
    %end
%end

% read dataset QC parameters if exist and override previous 
    % parameters file
    currentQCtest = mfilename;
    %ea_thresh_in_air = readDatasetParameter(sample_data.toolbox_input_file, currentQCtest, 'ea_thresh_in_air', ea_thresh_in_air);
    %if ea_thresh_in_air>20 % Reset value to default always to 20 counts
    %    ea_thresh_in_air=20;
    %end
% %% Input from user if no info is in the DDB
if ~isempty(noiseFloor)
    %noiseFloor = str2double(noiseFloor);
    ea_thresh_in_air = noiseFloor;
    
elseif isempty(noiseFloor)
    % read in filter parameters
    propFile  = fullfile('AutomaticQC', 'imosEchoIntensityInAirManualSetQC.txt');
    ea_thresh_in_air = str2double(readProperty('ea_thresh_in_air',   propFile));

    % read dataset QC parameters if exist and override previous 
    % parameters file
    %currentQCtest = mfilename;
    ea_thresh_in_air = readDatasetParameter(sample_data.toolbox_input_file, currentQCtest, 'ea_thresh_in_air', ea_thresh_in_air);
    if ea_thresh_in_air>20 % Reset value to default always to 20 counts
        ea_thresh_in_air=20;
    end
end

sizeCur = size(sample_data.variables{idUcur}.flags);

% Identify out-of-water period with 60 minute buffer
time_in_water = sample_data.time_deployment_start;
time_out_water = sample_data.time_deployment_end;
outofwater = false(sizeCur);
iOut = time <= time_in_water-(1/24);
iOut = iOut | time >= time_out_water+(1/24);
outofwater(iOut,:)=true;

% Verify if we have enough out-of-water data to compare with default value,
% if yes, ask the user if he wants to update ea_thresh % Aforest 16 February 2017
ea_p=(permute(ea,[2 3 1]));
   if beams ==4 %RDI 4 beams
        nam=['Echo intensity in air of ',sample_data.meta.instrument_model,' #',sample_data.meta.instrument_serial_no,' on ',sample_data.deployment_code];
   else %Nortek 3 beams
        nam=['Echo intensity in air of ',sample_data.meta.instrument_model,' #',sample_data.meta.instrument_serial_no,' on ',sample_data.deployment_code];
   end
if sum(sum(outofwater))==0
    ea_p=ea_p(:); %in counts
   
    % can be commented out as code requires additional toolbox licenses for nanmin
    ea_min=[sprintf('%.0f',nanmin(ea_p)),' counts'];
    ea_max=[sprintf('%.0f',nanmax(ea_p)),' counts'];
    ea_mean=[sprintf('%.0f',nanmean(ea_p)),' counts'];
    ea_median=[sprintf('%.0f',nanmedian(ea_p)),' counts'];
    ea_std=[sprintf('%.0f',nanstd(ea_p)),' counts'];
    %ea_p95=[sprintf('%.0f',(ea_p(round(0.95*length(ea_p))))),' counts'];
    
	%ea_p(ea_p<=0)=NaN; ea_p(isnan(ea_p))=[]; ea_p=sort(ea_p);% If any negative or NaN ea
	%ea_min=[sprintf('%.0f',min(ea_p)),' counts'];
	%ea_max=[sprintf('%.0f',max(ea_p)),' counts'];
	%ea_mean=[sprintf('%.0f',mean(ea_p)),' counts'];
	%ea_median=[sprintf('%.0f',median(ea_p)),' counts'];
	        
    str1=sprintf(['Info: imosEchoIntensityInAirManualSetQC cannot filter data based on echo intensity in air.\n\n***** Not enough out-of-water data available.']);
    str2=sprintf(['\n\nDo you want to use a threshold echo intensity flag based on the following echo intensity statistics when deployed?\n\nMin: ' ea_min '\nMax: ' ea_max '\nMean: ' ea_mean '\nMedian: ' ea_median '\nStandard deviation: ' ea_std]);
    str3=sprintf(['\n\nIf yes, please enter a given value. If not, click cancel.\n']);
    str4=[str1,str2,str3];
    xx = inputdlg(str4,nam,[1 120]);
    
    if ~isempty(xx)
    ea_thresh_in_air = str2num(xx{:}); % new ea_thresh_in_air 
    end
else
    ea_p=(ea_p(repmat(outofwater,[1,1,size(ea_p,3)]))); %in counts
      
    ea_min=[sprintf('%.0f',nanmin(ea_p)),' counts'];
    ea_max=[sprintf('%.0f',nanmax(ea_p)),' counts'];
    ea_mean=[sprintf('%.0f',nanmean(ea_p)),' counts'];
    ea_median=[sprintf('%.0f',nanmedian(ea_p)),' counts'];
    ea_std=[sprintf('%.0f',nanstd(ea_p)),' counts'];
   
   %ea_p(ea_p<=0)=NaN; ea_p(isnan(ea_p))=[]; ea_p=sort(ea_p);% If any negative or NaN ea
   %ea_min=[sprintf('%.0f',min(ea_p)),' counts'];
   %ea_max=[sprintf('%.0f',max(ea_p)),' counts'];
   %ea_mean=[sprintf('%.0f',mean(ea_p)),' counts'];
   %ea_median=[sprintf('%.0f',median(ea_p)),' counts'];
   %ea_std=[sprintf('%.0f',std(ea_p)),' counts'];
   %ea_p95=[sprintf('%.0f',ea_p(round(0.95*length(ea_p)))),' counts'];  %works only if there are no NANs, otherwise gives logical error
    
   str1=sprintf(['Info: imosEchoIntensityInAirManualSetQC requires your input\n\nThe following statistics for echo intensity in air have been found:\n\nMin: ' ea_min '\nMax: ' ea_max '\nMean: ' ea_mean '\nMedian: ' ea_median '\nStandard deviation: ' ea_std]);   
   str2=sprintf(['\n\nDo you want to use the database noted echo intensity  of ',int2str(ea_thresh_in_air), ' counts as criterion to flag velocity data?\n\nIf not, please enter a new value. If yes, click cancel.\n']); 
   str3=sprintf(['\nPlease note that using the mean echo in air is suggested, but using at least the minimum value is strongly advisable.\n']);
   str4=[str1,str2,str3];
   xx = inputdlg(str4,nam,[1 120],{int2str(ea_thresh_in_air)});
   
   if ~isempty(xx)
   ea_thresh_in_air = str2num(xx{:}); %ea_thresh_in_air  
   end
end

paramsLog = ['ea_thresh_in_air=' num2str(ea_thresh_in_air)];

%% same flags are given to any variable
flags = ones(sizeCur, 'int8')*rawFlag;

%% Run QC
% this test looks at each of the beam echo intensity and compares with
% ea_thresh_in_air value. The test passes if at least 3 beams (RDI) or 2
% beams (Nortek) provide echo intensity values above the ea_thresh_in_air
% value. Default is only 20 counts.
ea_p=permute(ea,[2 3 1]);
mask=ea_p>ea_thresh_in_air;
mask=sum(mask,3);
iFail=mask<beams-1;
iPass = ~iFail;

%% Run QC filter (iFail) on velocity data
flags(iFail) = badFlag;
flags(iPass) = goodFlag;

sample_data.variables{idUcur}.flags = flags;
sample_data.variables{idVcur}.flags = flags;
sample_data.variables{idWcur}.flags = flags;

varChecked = {sample_data.variables{idUcur}.name, ...
    sample_data.variables{idVcur}.name, ...
    sample_data.variables{idWcur}.name};

if idCdir
    sample_data.variables{idCdir}.flags = flags;
    varChecked = [varChecked, {sample_data.variables{idCdir}.name}];
end

if idCspd
    sample_data.variables{idCspd}.flags = flags;
    varChecked = [varChecked, {sample_data.variables{idCspd}.name}];
end

%% write/update dataset QC parameters
writeDatasetParameter(sample_data.toolbox_input_file, currentQCtest, 'ea_thresh_in_air', ea_thresh_in_air);

end
