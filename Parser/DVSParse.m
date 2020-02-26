function sample_data = DVSParse( filename, tMode )
%DVSPARSE Parses a raw (binary) data file from a Teledyne RD Workhorse DVS
% ADCP.
%
% This function uses the readDVSEnsembles function to read in a set
% of ensembles from a raw binary PD0 Workhorse ADCP file. It parses the 
% ensembles, and extracts and returns the following:
%
%   - time
%   - temperature (at each time)
%   - pressure (at each time, if present)
%   - salinity (at each time, if present)
%   - water speed (at each time and distance)
%   - water direction (at each time and distance)
%   - Acoustic backscatter intensity (at each time and distance, a separate 
%     variable for each beam)
%
% The conversion from the ADCP velocity values currently assumes that the 
% ADCP is using earth coordinates (see section 13.4 'Velocity Data Format' 
% of the Workhorse H-ADCP Operation Manual).
% 
% Inputs:
%   filename    - raw binary data file retrieved from a Workhorse.
%   tMode       - Toolbox data type mode.
%
% Outputs:
%   sample_data - sample_data struct containing the data retrieved from the
%                 input file.
%
% Author:       Shawn Meredyk <shawn.meredyk@as.ulaval.ca>
% Contributors: Leeying Wu <Wu.Leeying@saugov.sa.gov.au>
%               Bradley Morris <b.morris@unsw.edu.au>
%               Charles James May 2010 <charles.james@sa.gov.au>
%               Guillaume Galibert <guillaume.galibert@utas.edu.au>
%               Paul McCarthy <paul.mccarthy@csiro.au>
%         		Alexandre Forest <alexandre.forest@as.ulaval.ca>
%
% Copyright (c) 2017, Amundsen Science & ArcticNet
% http://www.amundsen.ulaval.ca/
% http://www.arcticnet.ulaval.ca/
% All rights reserved.
%
% Copyright (c) 2017, Australian Ocean Data Network (AODN) and Integrated 
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

  filename = filename{1};

  % we first look if the file has been processed to extract current data
  [filePath, fileRadName, ~] = fileparts(filename);
  
  currentFile   = fullfile(filePath, [fileRadName '.PD0']);
    
  
%  filename=which(filename); %prevent error with RDI files that are not found
%  if isempty(filename) 
%      filename=uigetfile([filePath,'\*.*'],['Select the file: ', fileRadName,'.000']); 
%      filename=[filePath,'\',filename];
%  end 
  
	ensembles = readDVSEnsembles( filename );
  
  if isempty(ensembles), error(['no ensembles found in file ' filename]); end
  
  %
  % retrieve metadata and data from struct
  %
  
  fixed = ensembles.fixedLeader;
  
  % metadata for this ensemble
  variable = ensembles.variableLeader;
  
  velocity = ensembles.velocity;
  
  backscatter1 = ensembles.echoIntensity.field1;
  backscatter2 = ensembles.echoIntensity.field2;
  backscatter3 = ensembles.echoIntensity.field3;
  backscatter4 = ensembles.echoIntensity.field4;
  
  correlation1 = ensembles.corrMag.field1;
  correlation2 = ensembles.corrMag.field2;
  correlation3 = ensembles.corrMag.field3;
  correlation4 = ensembles.corrMag.field4;
  
  percentGood1 = ensembles.percentGood.field1;
  percentGood2 = ensembles.percentGood.field2;
  percentGood3 = ensembles.percentGood.field3;
  percentGood4 = ensembles.percentGood.field4;
  clear ensembles;
  
  % we use these to set up variables and dimensions
  % we set a static value for these variables to the most frequent value found
  numBeams   = mode(fixed.numBeams);
  numCells   = mode(fixed.numCells);
  cellLength = mode(fixed.depthCellLength);
  cellStart  = mode(fixed.bin1Distance);
  
  % we can populate distance data now using cellLength and cellStart
  % ( / 100.0, as the ADCP gives the values in centimetres)
  cellStart  = cellStart  / 100.0;
  cellLength = cellLength / 100.0;
  
  % note this is actually distance between the ADCP's transducers and the
  % middle of each cell
  distance = (cellStart:  ...
      cellLength: ...
      cellStart + (numCells-1) * cellLength)';
  
  % rearrange the sample data
  instrument_firmware = strcat(num2str(fixed.cpuFirmwareVersion(1)), '.', num2str(fixed.cpuFirmwareRevision(1))); % we assume the first value is correct for the rest of the dataset
  %if str2double(instrument_firmware) > 8.35 % really meant for Workhorse units, as DVS units firmware is in the 40's anyway.
  %    time = datenum(...
  %        [variable.y2kCentury*100 + variable.y2kYear,...
  %        variable.y2kMonth,...
  %        variable.y2kDay,...
  %        variable.y2kHour,...
  %        variable.y2kMinute,...
  %        variable.y2kSecond + variable.y2kHundredth/100.0]);
  %else
      % looks like DVS Y2K compliant RTC time was not implemented
      century = 2000; 
      time = datenum(...
          century + variable.rtcYear,...
          variable.rtcMonth,...
          variable.rtcDay,...
          variable.rtcHour,...
          variable.rtcMinute,...
          variable.rtcSecond + variable.rtcHundredths/100.0);
    
  timePerPing = fixed.tppMinutes*60 + fixed.tppSeconds + fixed.tppHundredths/100;
  timePerEnsemble = fixed.pingsPerEnsemble .* timePerPing;
%   % shift the timestamp to the middle of the burst
%   time = time + (timePerEnsemble / (3600 * 24))/2;
%
 % try to guess model information
  adcpFreqs = str2num(fixed.systemConfiguration(:, 6:8)); % str2num is actually more relevant than str2double here
  adcpFreq = mode(adcpFreqs); % hopefully the most frequent value reflects the frequency when deployed
  switch adcpFreq
      case 0
          adcpFreq = 75;
          model = 'LongRanger';
		  xmtVolt = 2092719;
          
      case 1
          adcpFreq = 150;
          model = 'QuarterMaster';
          xmtVolt = 592157;
		  
      case 10
          adcpFreq = 300;
          model = 'Sentinel or Monitor';
          xmtVolt = 591257;
		  
      case 11
          adcpFreq = 600;
          model = 'Sentinel or Monitor';
          xmtVolt = 380667;
		  
      case 100
          adcpFreq = 1200;
          model = 'Sentinel or Monitor';
          xmtVolt = 253765;
		  
      otherwise
          adcpFreq = 2400;
          model = 'DVS';
          xmtVolt = 253765;
  end  
  %
  % auxillary data
  %
  temperature = variable.temperature;
 % pressure    = variable.pressure; %no pressure sensor in a DVS unit
%  salinity    = variable.salinity; %should let Toolbox calculate this
  pitch       = variable.pitch;
  roll        = variable.roll;
  heading     = variable.heading;
  voltageCnts = variable.adcChannel1;
  clear variable;
  
  %
  % calculate velocity (speed and direction)
  % currently assuming earth coordinate transform
  %
  
  veast = velocity.velocity1;
  vnrth = velocity.velocity2;
  wvel  = velocity.velocity3;
  evel  = velocity.velocity4;
  clear velocity;
  
  % set all bad values to NaN.
  vnrth(vnrth == -32768) = NaN;
  veast(veast == -32768) = NaN;
  wvel (wvel  == -32768) = NaN;
  evel (evel  == -32768) = NaN;
  
  %
  % temperature / 100.0  (0.01 deg   -> deg)
  % pressure    / 1000.0 (decapascal -> decibar)
  % vnrth       / 1000.0 (mm/s       -> m/s)
  % veast       / 1000.0 (mm/s       -> m/s)
  % wvel        / 1000.0 (mm/s       -> m/s)
  % evel        / 1000.0 (mm/s       -> m/s)
  % pitch       / 100.0  (0.01 deg   -> deg)
  % roll        / 100.0  (0.01 deg   -> deg)
  % heading     / 100.0  (0.01 deg   -> deg)
  % voltage 	/1000000.0 (counts   -> volt)
   % 
  %xmt voltage conversion for diagnostics
  % converting xmt voltage counts to volts , these are rough values
  
  % Voltage data tends to have many NaNs, therefore set all NaN to the previous value before it. 
  % Dimensions
[~,numCol] = size(voltageCnts);

% First, datai is copy of voltage data
datai = voltageCnts;

% For each column
for c = 1:numCol
    % Find first non-NaN row
    indxFirst = find(~isnan(voltageCnts(:,c)),1,'first');
    %if whole column is NaN
    %if( ~isempty(indxFirst) )
    % Find all NaN rows
    indxNaN = find(isnan(voltageCnts(:,c)));
    % Find NaN rows beyond first non-NaN
    indx = indxNaN(indxNaN > indxFirst);
    % For each of these, copy previous value
    for r = (indx(:))'
        datai(r,c) = datai(r-1,c);
    end
end    
  
    voltage	    = (datai*xmtVolt) /1000000; %  xmt voltage conversion , 
	%from p.136 of Workhorse Commands and Output Data Format PDF (RDI website - March 2016)
    
    clear datai voltageCnts numRow numCol indxFirst indxNaN indx r c;
    
  temperature  = temperature  / 100.0; 
  %pressure     = pressure     / 1000.0; %no pressure sensor in a DVS unit
  vnrth        = vnrth        / 1000.0;
  veast        = veast        / 1000.0;
  wvel         = wvel         / 1000.0;
  evel         = evel         / 1000.0;
  pitch        = pitch        / 100.0;
  roll         = roll         / 100.0;
  heading      = heading      / 100.0;
  
  % check for electrical/magnetic heading bias (usually magnetic declination)
  isMagBias = false;
  % we set a static value for this variable to the most frequent value found
  magDec = mode(fixed.headingBias)*0.01; % Scaling: LSD = 0.01degree; Range = -179.99 to 180.00degrees
  if magDec ~= 0
      isMagBias = true;
      magBiasComment = ['A compass correction of ' num2str(magDec) ...
          'degrees has been applied to the data by a technician using RDI''s software ' ...
          '(usually to account for magnetic declination).'];
  end
  
  speed = sqrt(vnrth.^2 + veast.^2);
  direction = getDirectionFromUV(veast, vnrth);

  
  % fill in the sample_data struct
  sample_data.toolbox_input_file        = filename;
  sample_data.meta.featureType          = ''; % strictly this dataset cannot be described as timeSeriesProfile since it also includes timeSeries data like TEMP
  sample_data.meta.fixedLeader          = fixed;
  sample_data.meta.binSize              = mode(fixed.depthCellLength)/100; % we set a static value for this variable to the most frequent value found
  sample_data.meta.instrument_make      = 'Teledyne RDI';
  
 
% this requires that the filename be properly annotated 'Make_Model_Serial' 
%extracting filename parts because DVS firmware doesn't record instrument
%Serial Number
unitInfo = textscan(fileRadName, '%s %s %s', 'Delimiter', '_');
serial = char(unitInfo{3});
  
  sample_data.meta.instrument_model             = [model ' ADCP'];
  sample_data.meta.instrument_serial_no         = serial;
  sample_data.meta.instrument_sample_interval   = median(diff(time*24*3600));
  sample_data.meta.instrument_average_interval  = mode(timePerEnsemble);
  sample_data.meta.instrument_firmware          = instrument_firmware;
  sample_data.meta.beam_angle              	 	= 45;  % DVS_OM_Jan09.pdf from the RDI support page, winADCP says 15 deg....
clear unitInfo serial;  % clear variables to avoid nomenclature conflicts, if any.
  
  % add dimensions with their data mapped
  adcpOrientations = str2num(fixed.systemConfiguration(:, 1)); % str2num is actually more relevant than str2double here
  adcpOrientation = mode(adcpOrientations); % hopefully the most frequent value reflects the orientation when deployed
  height = distance;
  if adcpOrientation == 0
      % case of a downward looking ADCP -> negative values
      height = -height;
      distance = -distance;
  end
  iBadOriented = adcpOrientations ~= adcpOrientation; % we'll only keep velocity data collected when ADCP is oriented as expected
  vnrth(iBadOriented, :) = NaN;
  veast(iBadOriented, :) = NaN;
  wvel(iBadOriented, :) = NaN;
  evel(iBadOriented, :) = NaN;
  speed(iBadOriented, :) = NaN;
  direction(iBadOriented, :) = NaN;
  backscatter1(iBadOriented, :) = NaN;
  backscatter2(iBadOriented, :) = NaN;
  backscatter3(iBadOriented, :) = NaN;
  backscatter4(iBadOriented, :) = NaN;
  correlation1(iBadOriented, :) = NaN;
  correlation2(iBadOriented, :) = NaN;
  correlation3(iBadOriented, :) = NaN;
  correlation4(iBadOriented, :) = NaN;
  percentGood1(iBadOriented, :) = NaN;
  percentGood2(iBadOriented, :) = NaN;
  percentGood3(iBadOriented, :) = NaN;
  percentGood4(iBadOriented, :) = NaN;
  dims = {
      'TIME',                   time,    ['Time stamp corresponds to the start of the measurement which lasts ' num2str(sample_data.meta.instrument_average_interval) ' seconds.']; ...
      'HEIGHT_ABOVE_SENSOR',    height,   'Data has been vertically bin-mapped using tilt information so that the cells have consistant heights above sensor in time.'; ...
      'DIST_ALONG_BEAMS',       distance, 'Data is not vertically bin-mapped (no tilt correction applied). Cells are lying parallel to the beams, at heights above sensor that vary with tilt.'
      };
  clear time height distance;
  
  nDims = size(dims, 1);
  sample_data.dimensions = cell(nDims, 1);
  for i=1:nDims
      sample_data.dimensions{i}.name         = dims{i, 1};
      sample_data.dimensions{i}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(dims{i, 1}, 'type')));
      sample_data.dimensions{i}.data         = sample_data.dimensions{i}.typeCastFunc(dims{i, 2});
      sample_data.dimensions{i}.comment      = dims{i, 3};
  end
  clear dims;
  
  % add information about the middle of the measurement period
  sample_data.dimensions{1}.seconds_to_middle_of_measurement = sample_data.meta.instrument_average_interval/2;
  
  % add variables with their dimensions and data mapped
  if isMagBias
      magExt = '';
  else
      magExt = '_MAG';
  end
  
  vars = {
      'TIMESERIES',         [],     1; ...
      'LATITUDE',           [],     NaN; ...
      'LONGITUDE',          [],     NaN; ...
      'NOMINAL_DEPTH',      [],     NaN; ...
      ['VCUR' magExt],      [1 2],  vnrth; ...
      ['UCUR' magExt],      [1 2],  veast; ...
      'WCUR',               [1 2],  wvel; ...
      'ECUR',               [1 2],  evel; ...
      'CSPD',               [1 2],  speed; ...
      ['CDIR' magExt],      [1 2],  direction; ...
      'ABSIC1',              [1 3],  backscatter1; ...
      'ABSIC2',              [1 3],  backscatter2; ...
      'ABSIC3',              [1 3],  backscatter3; ...
      'ABSIC4',              [1 3],  backscatter4; ...
      %'PRES_REL',            1,      pressure; ... %no pressure sensor for DVS units
      %'PSAL',               1,      salinity; ... % manual entry, could be incorrect, thus not used.
      'CMAG1',              [1 3],  correlation1; ...
      'CMAG2',              [1 3],  correlation2; ...
      'CMAG3',              [1 3],  correlation3; ...
      'CMAG4',              [1 3],  correlation4; ...
      'PERG1',              [1 2],  percentGood1; ...
      'PERG2',              [1 2],  percentGood2; ...
      'PERG3',              [1 2],  percentGood3; ...
      'PERG4',              [1 2],  percentGood4; ...
      'PITCH',              1,      pitch; ...
      'ROLL',               1,      roll; ...
	  'TEMP',               1,      temperature; ...
	  'VOLT',				1,		voltage; ... % added for equipment diagnostics 
      ['HEADING' magExt],   1,      heading
      };
  
  clear vnrth veast wvel evel speed direction backscatter1 ...
      backscatter2 backscatter3 backscatter4 temperature ... % remove pressure from the list, as DVS has no pressure sensor
      correlation1 correlation2 correlation3 correlation4 ... % salinity removed from the list, let the toolbox calculate this
      percentGood1 percentGood2 percentGood3 percentGood4 pitch roll ...
      heading voltage; % added voltage to the list
  
  nVars = size(vars, 1);
  sample_data.variables = cell(nVars, 1);
  for i=1:nVars
      sample_data.variables{i}.name         = vars{i, 1};
      sample_data.variables{i}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(vars{i, 1}, 'type')));
      sample_data.variables{i}.dimensions   = vars{i, 2};
      
      % we don't want coordinates attribute for LATITUDE, LONGITUDE and NOMINAL_DEPTH
      if ~isempty(sample_data.variables{i}.dimensions)
          switch sample_data.variables{i}.dimensions(end)
              case 1
                  sample_data.variables{i}.coordinates = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
              case 2
                  sample_data.variables{i}.coordinates = 'TIME LATITUDE LONGITUDE HEIGHT_ABOVE_SENSOR';
              case 3
                  sample_data.variables{i}.coordinates = 'TIME LATITUDE LONGITUDE DIST_ALONG_BEAMS';
          end
      end
      
      sample_data.variables{i}.data         = sample_data.variables{i}.typeCastFunc(vars{i, 3});
      
      % no pressure sensor, so below text is not needed.
      %if strcmpi(vars{i, 1}, 'PRES_REL')
      %    sample_data.variables{i}.applied_offset = sample_data.variables{i}.typeCastFunc(-gsw_P0/10^4); % (gsw_P0/10^4 = 10.1325 dbar)
      %end
      
      if any(strcmpi(vars{i, 1}, {'VCUR', 'UCUR', 'CDIR', 'HEADING'}))
          sample_data.variables{i}.compass_correction_applied = magDec;
          sample_data.variables{i}.comment = magBiasComment;
      end
  end
  clear vars;
  
  % remove auxillary data if the sensors 
  % were not installed on the instrument
  hasPres    = mode(str2num(fixed.sensorsAvailable(:, 3))); % str2num is actually more relevant than str2double here,also DVS units do not have pressure sensors
  hasHeading = mode(str2num(fixed.sensorsAvailable(:, 4)));
  hasPitch   = mode(str2num(fixed.sensorsAvailable(:, 5)));
  hasRoll    = mode(str2num(fixed.sensorsAvailable(:, 6)));
  hasPsal    = mode(str2num(fixed.sensorsAvailable(:, 7)));
  hasTemp    = mode(str2num(fixed.sensorsAvailable(:, 8)));

  % indices of variables to remove
  remove = [];
  
  if ~hasPres,    remove(end+1) = getVar(sample_data.variables, 'PRES_REL');end
  if ~hasHeading, remove(end+1) = getVar(sample_data.variables, 'HEADING'); end
  if ~hasPitch,   remove(end+1) = getVar(sample_data.variables, 'PITCH');   end
  if ~hasRoll,    remove(end+1) = getVar(sample_data.variables, 'ROLL');    end
  if ~hasPsal,    remove(end+1) = getVar(sample_data.variables, 'PSAL');    end
  if ~hasTemp,    remove(end+1) = getVar(sample_data.variables, 'TEMP');    end
  
  % also remove empty backscatter and correlation data in case of ADCP with
  % less than 4 beams
  for k = 4:-1:numBeams+1
      kStr = num2str(k);
      remove(end+1) = getVar(sample_data.variables, ['ABSIC' kStr]);
      remove(end+1) = getVar(sample_data.variables, ['CMAG' kStr]);
      remove(end+1) = getVar(sample_data.variables, ['PERG' kStr]);
  end 
 
end % end of main function

function direction = getDirectionFromUV(uvel, vvel)
    % direction is in degrees clockwise from north
    direction = atan(abs(uvel ./ vvel)) .* (180 / pi);
    
    % !!! if vvel == 0 we get NaN !!!
    direction(vvel == 0) = 90;
    
    se = vvel <  0 & uvel >= 0;
    sw = vvel <  0 & uvel <  0;
    nw = vvel >= 0 & uvel <  0;
    
    direction(se) = 180 - direction(se);
    direction(sw) = 180 + direction(sw);
    direction(nw) = 360 - direction(nw);
end