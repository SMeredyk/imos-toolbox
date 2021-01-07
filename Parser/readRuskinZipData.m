function sample_data = readRSKfile( filename, mode )
%% Function Description
% readRSKfile.m Parses a data file retrieved from an RBR logger 
% that was exported from the new (fall 2020) ruskin v2.12.
% this script will extract the data needed to process an .rsk file
% without having to manually export the Ruskin Zip folder.
%
% Inputs:
%   filename    - Cell array containing the name of the file to parse.
%   mode        - Toolbox data type mode ('profile' or 'timeSeries').
%
% Outputs:
%   sample_data - Struct containing imported sample data.
%
% Author : 		Shawn Meredyk <shawn.meredyk@arcticnet.ulaval.ca>
% Contributor : Guillaume Galibert <guillaume.galibert@utas.edu.au>
%				
%
% Copyright (c) 2020, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS) and Amundsen Science (AS).
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
% End of Function Description

  narginchk(2,2);
  
  if ~ischar(filename)  
    error('filename must be a string'); 
  end
  
  % open the file, and read in the header and data
  try 
    
    fid    = fopen(filename, 'rt');
    header = readHeader(fid);
    data   = readData(fid, header);
    fclose(fid);
  
  catch e
    if fid ~= -1, fclose(fid); end
    rethrow(e);
  end
  
  % copy all of the information over to the sample data struct
  sample_data = struct;

  sample_data.toolbox_input_file                = filename;
  sample_data.meta.instrument_make              = header.make;
  sample_data.meta.instrument_model             = header.model;
  sample_data.meta.instrument_firmware          = header.firmware;
  sample_data.meta.instrument_serial_no         = header.serial;
  sample_data.meta.instrument_sample_interval   = median(diff(data.time*24*3600));
  sample_data.meta.featureType                  = mode;
  
  sample_data.dimensions = {};  
  sample_data.variables  = {};
  
    %% Case Profile
  switch mode
      case 'profile'
     
          % dimensions creation
          iVarPRES = NaN;
          iVarDEPTH = NaN;
          isZ = false;
          vars = fieldnames(data);
          nVars = length(vars);
          for k = 1:nVars
              if strcmpi('DEPTH', vars{k})
                  iVarDEPTH = k;
                  isZ = true;
                  break;
              end
              if strcmpi('PRES', vars{k})
                  iVarPRES = k;
                  isZ = true;
              end
              if ~isnan(iVarDEPTH) && ~isnan(iVarPRES), break; end
          end
          
          if ~isZ
              error('There is no pressure or depth information in this file to use it in profile mode');
          end
          
          depthComment = '';
          if ~isnan(iVarDEPTH)
              iVarZ = iVarDEPTH;
              depthData = data.(vars{iVarDEPTH});
          else
              iVarZ = iVarPRES;
              depthData = data.(vars{iVarPRES} - gsw_P0/10^4);
              presComment = ['abolute '...
                  'pressure measurements to which a nominal '...
                  'value for atmospheric pressure (10.1325 dbar) '...
                  'has been substracted'];
              depthComment  = ['Depth computed from '...
                  presComment ', assuming 1dbar ~= 1m.'];
          end
          
          % let's distinguish descending/ascending parts of the profile
          nData = length(data.(vars{iVarZ}));
          zMax = max(data.(vars{iVarZ}));
          posZMax = find(data.(vars{iVarZ}) == zMax, 1, 'last'); % in case there are many times the max value
          iD = [true(posZMax, 1); false(nData-posZMax, 1)];
          
          nD = sum(iD);
          nA = sum(~iD);
          MAXZ = max(nD, nA);
          
          dNaN = nan(MAXZ-nD, 1);
          aNaN = nan(MAXZ-nA, 1);
          
          if nA == 0
              sample_data.dimensions{1}.name            = 'DEPTH';
              sample_data.dimensions{1}.typeCastFunc    = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{1}.name, 'type')));
              sample_data.dimensions{1}.data            = sample_data.dimensions{1}.typeCastFunc(depthData);
              sample_data.dimensions{1}.comment         = depthComment;
              sample_data.dimensions{1}.axis            = 'Z';
              
              sample_data.variables{end+1}.name         = 'PROFILE';
              sample_data.variables{end}.typeCastFunc   = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
              sample_data.variables{end}.data           = sample_data.variables{end}.typeCastFunc(1);
              sample_data.variables{end}.dimensions     = [];
          else
              sample_data.dimensions{1}.name            = 'MAXZ';
              sample_data.dimensions{1}.typeCastFunc    = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{1}.name, 'type')));
              sample_data.dimensions{1}.data            = sample_data.dimensions{1}.typeCastFunc(1:1:MAXZ);
              
              sample_data.dimensions{2}.name            = 'PROFILE';
              sample_data.dimensions{2}.typeCastFunc    = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{2}.name, 'type')));
              sample_data.dimensions{2}.data            = sample_data.dimensions{2}.typeCastFunc([1, 2]);
              
              disp(['Warning : ' sample_data.toolbox_input_file ...
                  ' is not IMOS CTD profile compliant. See ' ...
                  'http://help.aodn.org.au/help/sites/help.aodn.org.au/' ...
                  'files/ANMN%20CTD%20Processing%20Procedures.pdf']);
          end
          
          % Add TIME, DIRECTION and POSITION infos
          descendingTime = data.time(iD);
          descendingTime = descendingTime(1);
          
          if nA == 0
              ascendingTime = [];
              dimensions = [];
          else
              ascendingTime = data.time(~iD);
              ascendingTime = ascendingTime(1);
              dimensions = 2;
          end
          
          sample_data.variables{end+1}.dimensions   = dimensions;
          sample_data.variables{end}.name         = 'TIME';
          sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
          sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc([descendingTime, ascendingTime]);
          sample_data.variables{end}.comment      = 'First value over profile measurement';
          
          sample_data.variables{end+1}.dimensions = dimensions;
          sample_data.variables{end}.name         = 'DIRECTION';
          sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
          if nA == 0
              sample_data.variables{end}.data     = {'D'};
          else
              sample_data.variables{end}.data     = {'D', 'A'};
          end
          
          sample_data.variables{end+1}.dimensions = dimensions;
          sample_data.variables{end}.name         = 'LATITUDE';
          sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
          if nA == 0
              sample_data.variables{end}.data     = sample_data.variables{end}.typeCastFunc(NaN);
          else
              sample_data.variables{end}.data     = sample_data.variables{end}.typeCastFunc([NaN, NaN]);
          end
          
          sample_data.variables{end+1}.dimensions = dimensions;
          sample_data.variables{end}.name         = 'LONGITUDE';
          sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
          if nA == 0
              sample_data.variables{end}.data     = sample_data.variables{end}.typeCastFunc(NaN);
          else
              sample_data.variables{end}.data     = sample_data.variables{end}.typeCastFunc([NaN, NaN]);
          end
          
          sample_data.variables{end+1}.dimensions = dimensions;
          sample_data.variables{end}.name         = 'BOT_DEPTH';
          sample_data.variables{end}.comment      = 'Bottom depth measured by ship-based acoustic sounder at time of CTD cast.';
          sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
          if nA == 0
              sample_data.variables{end}.data     = sample_data.variables{end}.typeCastFunc(NaN);
          else
              sample_data.variables{end}.data     = sample_data.variables{end}.typeCastFunc([NaN, NaN]);
          end
          
          % Manually add variable DEPTH if multiprofile and doesn't exit
          % yet
          if isnan(iVarDEPTH) && (nA ~= 0)
              sample_data.variables{end+1}.dimensions = [1 2];
              
              sample_data.variables{end}.name         = 'DEPTH';
              sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
              
              % we need to padd data with NaNs so that we fill MAXZ
              % dimension
              sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc([[depthData(iD); dNaN], [depthData(~iD); aNaN]]);
              
              sample_data.variables{end}.comment      = depthComment;
              sample_data.variables{end}.axis         = 'Z';
          end
          
          % scan through the list of parameters that were read
          % from the file, and create a variable for each
          for k = 1:nVars
              % we skip TIME and DEPTH
              if strcmpi('TIME', vars{k}), continue; end
              if strcmpi('DEPTH', vars{k}) && (nA == 0), continue; end
              
              name = '';
              comment.(vars{k}) = '';
              switch vars{k}
                  
                  %Conductivity (mS/cm) = 10-1*(S/m)
                  case {'Cond', 'cond00', 'cond05', 'cond06'},
                      name = 'CNDC';
                      data.(fields{k}) = data.(fields{k})/10;
                      
                      %Temperature (Celsius degree)
                  case {'Temp', 'temp00', 'temp09'}, name = 'TEMP';
                      
                      %Pressure (dBar)
                  case {'Pres', 'pres19', 'pres24'}, name = 'PRES';
                      
                      %Pressure_Relative (dBar)
                  case 'pres08', name = 'PRES_REL';
                      
                      %Fluorometry-chlorophyl (ug/l) = (mg.m-3)
                  case {'FlC', 'fluo01'},
                      name = 'CPHL';
                      comment.(fields{k}) = ['Artificial chlorophyll data computed from ' ...
                          'fluorometry sensor raw counts measurements. Originally ' ...
                          'expressed in ug/l, 1l = 0.001m3 was assumed.'];
                      
                      %Turbidity (NTU)
                  case {'Turb', 'Turb-a','turb00'}, name = 'TURB';
                      comment.(fields{k}) = 'NTU - Auto Ranging Seapoint SCF';
                      
                      %Rinko temperature (Celsius degree)
                  case 'R_Temp',
                      name = '';
                      comment.(fields{k}) = 'Corrected temperature.';
                      
                      %Rinko dissolved O2 (%)
                  case 'R_D_O2', name = 'DOXS';
                      
                      %Depth (m)
                  % case {'Depth', 'dpth01'}, name = 'DEPTH'; % letting toolbox calculate this.
                      
                      %Salinity (PSU)
                  %case 'sal_00', name = 'PSAL'; % removed 'Salin' due to possible bad salinty data from Ruskin
                      
                      %Specific conductivity (uS/cm) = 10-4 * (S/m)
                  case {'SpecCond', 'specC', 'scon00'},
                      name = 'SPEC_CNDC';
                      data.(fields{k}) = data.(fields{k})/10000;
                      
                      %Density anomaly (n/a)
                  %case {'DensAnom', 'density', 'dden00'}, name = '';
                      
                      %Speed of sound (m/s)
                  %case 'sos_00', name = 'SSPD'; % removed 'SosUN' due to possible bad salinty data from Ruskin
                      
                      %Rinko dissolved O2 concentration (mg/l) => (umol/l)
                  case 'rdO2C',
                      name = 'DOX1';
                      comment.(fields{k}) = ['Originally expressed in mg/l, ' ...
                          'O2 density = 1.429kg/m3 and 1ml/l = 44.660umol/l were assumed.'];
                      data.(fields{k}) = data.(fields{k}) * 44.660/1.429; % O2 density = 1.429 kg/m3
                      
                      % Oxyguard dissolved O2 (%)
                  case 'D_O2', name = 'DOXS';
                      
                      % Oxyguard dissolved O2 concentration (ml/l) => (umol/l)
                  case {'dO2C','xdO2C'},
                      name = 'DOX1';
                      comment.(fields{k}) = ['Originally expressed in ml/l, ' ...
                          '1ml/l = 44.660umol/l was assumed.'];
                      data.(fields{k}) = data.(fields{k}) * 44.660;
              end
              
              if ~isempty(name)
                  sample_data.variables{end+1}.dimensions = [1 dimensions];
                  
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  if nA == 0
                      sample_data.variables{end  }.data   = sample_data.variables{end}.typeCastFunc(data.(vars{k})(iD));
                  else
                      % we need to padd data with NaNs so that we fill MAXZ
                      % dimension
                      sample_data.variables{end  }.data   = sample_data.variables{end}.typeCastFunc([[data.(vars{k})(iD); dNaN], [data.(vars{k})(~iD); aNaN]]);
                  end
                  sample_data.variables{end}.comment    = comment.(vars{k});
                  
                  if ~any(strcmpi(vars{k}, {'TIME', 'DEPTH'}))
                      sample_data.variables{end  }.coordinates = 'TIME LATITUDE LONGITUDE DEPTH';
                  end
              end   
          end
          
          % Let's add DOX1/DOX2 if PSAL/CNDC, TEMP and DOXS are present and DOX1 not
          % already present
          doxs = getVar(sample_data.variables, 'DOXS');
          dox1 = getVar(sample_data.variables, 'DOX1');
          if doxs ~= 0 && dox1 == 0
              doxs = sample_data.variables{doxs};
              name = 'DOX1';
              
              % to perform this conversion, we need temperature,
              % and salinity/conductivity+pressure data to be present
              temp = getVar(sample_data.variables, 'TEMP');
              psal = getVar(sample_data.variables, 'PSAL');
              cndc = getVar(sample_data.variables, 'CNDC');
              pres = getVar(sample_data.variables, 'PRES');
              
              % if any of this data isn't present,
              % we can't perform the conversion
              if temp ~= 0 && (psal ~= 0 || (cndc ~= 0 && pres ~= 0))
                  temp = sample_data.variables{temp};
                  if psal ~= 0
                      psal = sample_data.variables{psal};
                  else
                      cndc = sample_data.variables{cndc};
                      pres = sample_data.variables{pres};
                      % conductivity is in S/m and gsw_C3515 in mS/cm
                      crat = 10*cndc.data ./ gsw_C3515;
                      
                      % we need to use relative pressure using gsw_P0 = 101325 Pa 
                      psal.data = gsw_SP_from_R(crat, temp.data, pres.data - gsw_P0/10^4);
                  end
                  
                  % O2 solubility (Garcia and Gordon, 1992-1993)
                  %
                  solubility = O2sol(psal.data, temp.data, 'ml/l');
                  
                  % O2 saturation to O2 concentration measured
                  % O2 saturation (per cent) = 100* [O2/O2sol]
                  %
                  % that is to say : O2 = O2sol * O2sat / 100
                  data = solubility .* doxs.data / 100;
                  
                  % conversion from ml/l to umol/l
                  data = data * 44.660;
                  comment = ['Originally expressed in % of saturation, using Garcia '...
                      'and Gordon equations (1992-1993) and ml/l coefficients, assuming 1ml/l = 44.660umol/l.'];
                  
                  sample_data.variables{end+1}.dimensions = [1 dimensions];
                  sample_data.variables{end}.comment      = comment;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  if nA == 0
                      sample_data.variables{end}.data   = sample_data.variables{end}.typeCastFunc(data(iD));
                  else
                      % we need to padd data with NaNs so that we fill MAXZ
                      % dimension
                      sample_data.variables{end}.data   = sample_data.variables{end}.typeCastFunc([[data(iD); dNaN], [data(~iD); aNaN]]);
                  end
                  sample_data.variables{end}.coordinates = 'TIME LATITUDE LONGITUDE DEPTH';
                  
                  % Let's add DOX2
                  name = 'DOX2';
                  
                  % O2 solubility (Garcia and Gordon, 1992-1993)
                  %
                  solubility = O2sol(psal.data, temp.data, 'umol/kg');
                  
                  % O2 saturation to O2 concentration measured
                  % O2 saturation (per cent) = 100* [O2/O2sol]
                  %
                  % that is to say : O2 = O2sol * O2sat / 100
                  data = solubility .* doxs.data / 100;
                  comment = ['Originally expressed in % of saturation, using Garcia '...
                      'and Gordon equations (1992-1993) and umol/kg coefficients.'];
                  
                  sample_data.variables{end+1}.dimensions = [1 dimensions];
                  sample_data.variables{end}.comment      = comment;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  if nA == 0
                      sample_data.variables{end}.data   = sample_data.variables{end}.typeCastFunc(data(iD));
                  else
                      % we need to padd data with NaNs so that we fill MAXZ
                      % dimension
                      sample_data.variables{end}.data   = sample_data.variables{end}.typeCastFunc([[data(iD); dNaN], [data(~iD); aNaN]]);
                  end
                  sample_data.variables{end}.coordinates = 'TIME LATITUDE LONGITUDE DEPTH';
              end
          end
          
          % Let's add a new parameter if DOX1, PSAL/CNDC, TEMP and PRES are
          % present and DOX2 not already present
          dox1 = getVar(sample_data.variables, 'DOX1');
          dox2 = getVar(sample_data.variables, 'DOX2');
          if dox1 ~= 0 && dox2 == 0
              dox1 = sample_data.variables{dox1};
              name = 'DOX2';
              
              % umol/l -> umol/kg
              %
              % to perform this conversion, we need to calculate the
              % density of sea water; for this, we need temperature,
              % salinity, and pressure data to be present
              temp = getVar(sample_data.variables, 'TEMP');
              pres = getVar(sample_data.variables, 'PRES');
              psal = getVar(sample_data.variables, 'PSAL');
              cndc = getVar(sample_data.variables, 'CNDC');
              
              % if any of this data isn't present,
              % we can't perform the conversion to umol/kg
              if temp ~= 0 && pres ~= 0 && (psal ~= 0 || cndc ~= 0)
                  temp = sample_data.variables{temp};
                  pres = sample_data.variables{pres};
                  if psal ~= 0
                      psal = sample_data.variables{psal};
                  else
                      cndc = sample_data.variables{cndc};
                      % conductivity is in S/m and gsw_C3515 in mS/cm
                      crat = 10*cndc.data ./ gsw_C3515;
                      
                      % we need to use relative pressure using gsw_P0 = 101325 Pa 
                      psal.data = gsw_SP_from_R(crat, temp.data, pres.data - gsw_P0/10^4);
                  end
                  
                  % calculate density from salinity, temperature and pressure
                  dens = sw_dens(psal.data, temp.data, pres.data - gsw_P0/10^4); % cannot use the GSW SeaWater library TEOS-10 as we don't know yet the position
                  
                  % umol/l -> umol/kg (dens in kg/m3 and 1 m3 = 1000 l)
                  data = dox1.data .* 1000.0 ./ dens;
                  comment = ['Originally expressed in mg/l, assuming O2 density = 1.429kg/m3, 1ml/l = 44.660umol/l '...
                      'and using density computed from Temperature, Salinity and Pressure '...
                      'with the CSIRO SeaWater library (EOS-80) v1.1.'];
                  
                  sample_data.variables{end+1}.dimensions = [1 dimensions];
                  sample_data.variables{end}.comment      = comment;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  if nA == 0
                      sample_data.variables{end}.data   = sample_data.variables{end}.typeCastFunc(data(iD));
                  else
                      % we need to padd data with NaNs so that we fill MAXZ
                      % dimension
                      sample_data.variables{end}.data   = sample_data.variables{end}.typeCastFunc([[data(iD); dNaN], [data(~iD); aNaN]]);
                  end
                  sample_data.variables{end}.coordinates = 'TIME LATITUDE LONGITUDE DEPTH';
              end
          end
      % End Case Profile    
      otherwise
      %% Case TimeSeries    
          sample_data.dimensions{1}.name            = 'TIME';
          sample_data.dimensions{1}.typeCastFunc    = str2func(netcdf3ToMatlabType(imosParameters(sample_data.dimensions{1}.name, 'type')));
          sample_data.dimensions{1}.data            = sample_data.dimensions{1}.typeCastFunc(data.time);
          
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
          data = rmfield(data, 'time');
          fields = fieldnames(data);
          coordinates = 'TIME LATITUDE LONGITUDE NOMINAL_DEPTH';
          
          for k = 1:length(fields)
              
              name = '';
              comment.(fields{k}) = '';
              switch fields{k}
                  
                  %Conductivity (mS/cm) = 10-1*(S/m)
                  case {'Cond', 'cond00', 'cond05', 'cond06'},
                      name = 'CNDC';
                      data.(fields{k}) = data.(fields{k})/10;
                      
                      %Temperature (Celsius degree)
                  case {'Temp', 'temp00', 'temp09'}, name = 'TEMP';
                      
                      %Pressure (dBar)
                  case {'Pres', 'pres19', 'pres24'}, name = 'PRES';
                      
                      %Pressure_Relative (dBar)
                  case 'pres08', name = 'PRES_REL';
                      
                      %Fluorometry-chlorophyl (ug/l) = (mg.m-3)
                  case {'FlC', 'fluo01'},
                      name = 'CPHL';
                      comment.(fields{k}) = ['Artificial chlorophyll data computed from ' ...
                          'fluorometry sensor raw counts measurements. Originally ' ...
                          'expressed in ug/l, 1l = 0.001m3 was assumed.'];
                      
                      %Turbidity (NTU)
                  case {'Turb', 'Turb-a','turb00'}, name = 'TURB';
                      comment.(fields{k}) = 'NTU - Auto Ranging Seapoint SCF';
                      
                      %Rinko temperature (Celsius degree)
                  case 'R_Temp',
                      name = '';
                      comment.(fields{k}) = 'Corrected temperature.';
                      
                      %Rinko dissolved O2 (%)
                  case 'R_D_O2', name = 'DOXS';
                      
                      %Depth (m)
                  %% case {'Depth', 'dpth01'}, name = 'DEPTH'; % letting the toolbox calculate Depth
                      
                      %Salinity (PSU)
                 %% case 'sal_00', name = 'PSAL'; % removed 'Salin' due to possible bad salinty data from Ruskin
                      
                      %Specific conductivity (uS/cm) = 10-4 * (S/m)
                  case {'SpecCond', 'specC', 'scon00'},
                      name = 'SPEC_CNDC';
                      data.(fields{k}) = data.(fields{k})/10000;
                      
                      %Density anomaly (n/a)
                 % case {'DensAnom', 'density', 'dden00'}, name = '';
                      
                      %Speed of sound (m/s)
                %%  case 'sos_00', name = 'SSPD'; % removed 'SosUN' due to possible bad salinty data from Ruskin
                      
                      %Rinko dissolved O2 concentration (mg/l) => (umol/l)
                  case 'rdO2C',
                      name = 'DOX1';
                      comment.(fields{k}) = ['Originally expressed in mg/l, ' ...
                          'O2 density = 1.429kg/m3 and 1ml/l = 44.660umol/l were assumed.'];
                      data.(fields{k}) = data.(fields{k}) * 44.660/1.429; % O2 density = 1.429 kg/m3
                      
                      % Oxyguard dissolved O2 (%)
                  case 'D_O2', name = 'DOXS';
                      
                      % Oxyguard dissolved O2 concentration (ml/l) => (umol/l)
                  case {'dO2C','xdO2C'},
                      name = 'DOX1';
                      comment.(fields{k}) = ['Originally expressed in ml/l, ' ...
                          '1ml/l = 44.660umol/l was assumed.'];
                      data.(fields{k}) = data.(fields{k}) * 44.660;
              end
              
              if ~isempty(name)
                  % dimensions definition must stay in this order : T, Z, Y, X, others;
                  % to be CF compliant
                  sample_data.variables{end+1}.dimensions = 1;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc(data.(fields{k}));
                  sample_data.variables{end}.coordinates  = coordinates;
                  sample_data.variables{end}.comment      = comment.(fields{k});
              end
          end
          
          % Let's add DOX1/DOX2 if PSAL/CNDC, TEMP and DOXS are present and DOX1 not
          % already present
          doxs = getVar(sample_data.variables, 'DOXS');
          dox1 = getVar(sample_data.variables, 'DOX1');
          if doxs ~= 0 && dox1 == 0
              doxs = sample_data.variables{doxs};
              name = 'DOX1';
              
              % to perform this conversion, we need temperature,
              % and salinity/conductivity+pressure data to be present
              temp = getVar(sample_data.variables, 'TEMP');
              psal = getVar(sample_data.variables, 'PSAL');
              cndc = getVar(sample_data.variables, 'CNDC');
              pres = getVar(sample_data.variables, 'PRES');
              
              % if any of this data isn't present,
              % we can't perform the conversion
              if temp ~= 0 && (psal ~= 0 || (cndc ~= 0 && pres ~= 0))
                  temp = sample_data.variables{temp};
                  if psal ~= 0
                      psal = sample_data.variables{psal};
                  else
                      cndc = sample_data.variables{cndc};
                      pres = sample_data.variables{pres};
                      % conductivity is in S/m and gsw_C3515 in mS/cm
                      crat = 10*cndc.data ./ gsw_C3515;
                      
                      % we need to use relative pressure using gsw_P0 = 101325 Pa 
                      psal.data = gsw_SP_from_R(crat, temp.data, pres.data - gsw_P0/10^4);
                  end
                  
                  % O2 solubility (Garcia and Gordon, 1992-1993)
                  %
                  solubility = O2sol(psal.data, temp.data, 'ml/l');
                  
                  % O2 saturation to O2 concentration measured
                  % O2 saturation (per cent) = 100* [O2/O2sol]
                  %
                  % that is to say : O2 = O2sol * O2sat / 100
                  data = solubility .* doxs.data / 100;
                  
                  % conversion from ml/l to umol/l
                  data = data * 44.660;
                  comment = ['Originally expressed in % of saturation, using Garcia '...
                      'and Gordon equations (1992-1993) and ml/l coefficients, assuming 1ml/l = 44.660umol/l.'];
                  
                  sample_data.variables{end+1}.dimensions = 1;
                  sample_data.variables{end}.comment      = comment;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc(data);
                  sample_data.variables{end}.coordinates  = coordinates;
                  
                  % Let's add DOX2
                  name = 'DOX2';
                  
                  % O2 solubility (Garcia and Gordon, 1992-1993)
                  %
                  solubility = O2sol(psal.data, temp.data, 'umol/kg');
                  
                  % O2 saturation to O2 concentration measured
                  % O2 saturation (per cent) = 100* [O2/O2sol]
                  %
                  % that is to say : O2 = O2sol * O2sat / 100
                  data = solubility .* doxs.data / 100;
                  comment = ['Originally expressed in % of saturation, using Garcia '...
                      'and Gordon equations (1992-1993) and umol/kg coefficients.'];
                  
                  sample_data.variables{end+1}.dimensions = 1;
                  sample_data.variables{end}.comment      = comment;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc(data);
                  sample_data.variables{end}.coordinates  = coordinates;
              end
          end
          
          % Let's add a new parameter if DOX1, PSAL/CNDC, TEMP and PRES are
          % present and DOX2 not already present
          dox1 = getVar(sample_data.variables, 'DOX1');
          dox2 = getVar(sample_data.variables, 'DOX2');
          if dox1 ~= 0 && dox2 == 0
              dox1 = sample_data.variables{dox1};
              name = 'DOX2';
              
              % umol/l -> umol/kg
              %
              % to perform this conversion, we need to calculate the
              % density of sea water; for this, we need temperature,
              % salinity, and pressure data to be present
              temp = getVar(sample_data.variables, 'TEMP');
              pres = getVar(sample_data.variables, 'PRES');
              psal = getVar(sample_data.variables, 'PSAL');
              cndc = getVar(sample_data.variables, 'CNDC');
              
              % if any of this data isn't present,
              % we can't perform the conversion to umol/kg
              if temp ~= 0 && pres ~= 0 && (psal ~= 0 || cndc ~= 0)
                  temp = sample_data.variables{temp};
                  pres = sample_data.variables{pres};
                  if psal ~= 0
                      psal = sample_data.variables{psal};
                  else
                      cndc = sample_data.variables{cndc};
                      % conductivity is in S/m and gsw_C3515 in mS/cm
                      crat = 10*cndc.data ./ gsw_C3515;
                      
                      % we need to use relative pressure using gsw_P0 = 101325 Pa 
                      psal.data = gsw_SP_from_R(crat, temp.data, pres.data - gsw_P0/10^4);
                  end
                  
                  % calculate density from salinity, temperature and pressure
                  dens = sw_dens(psal.data, temp.data, pres.data - gsw_P0/10^4); % cannot use the GSW SeaWater library TEOS-10 as we don't know yet the position
                  
                  % umol/l -> umol/kg (dens in kg/m3 and 1 m3 = 1000 l)
                  data = dox1.data .* 1000.0 ./ dens;
                  comment = ['Originally expressed in mg/l, assuming O2 density = 1.429kg/m3, 1ml/l = 44.660umol/l '...
                      'and using density computed from Temperature, Salinity and Pressure '...
                      'with the CSIRO SeaWater library (EOS-80) v1.1.'];
                  
                  sample_data.variables{end+1}.dimensions = 1;
                  sample_data.variables{end}.comment      = comment;
                  sample_data.variables{end}.name         = name;
                  sample_data.variables{end}.typeCastFunc = str2func(netcdf3ToMatlabType(imosParameters(sample_data.variables{end}.name, 'type')));
                  sample_data.variables{end}.data         = sample_data.variables{end}.typeCastFunc(data);
                  sample_data.variables{end}.coordinates  = coordinates;
              end
          end
  end
end
% End Case TimeSeries
  
function header = readHeader(fid)
%READHEADER Reads the header section from the top of the file.

  header = struct;
  lines  = {};
  
  line = fgetl(fid);
  
  while ~contains(line, 'Date & Time')
    lines = [lines line];
    line  = fgetl(fid);
  end
  
  header.variables = strtrim(line);
  
  % use regexp to read in all the important header information
  exprs = {
    ['^Model=+' '(\S+)$']
    ['^Firmware=+' '(\S+)$']
    ['^Serial=+' '(\S+)$']
    ['^LoggingStartDate=+' '(\S+)$'] % this variable doesn`t exist in RBRConcerto files 
    ['^LoggingStartTime=+' '(\S+\s?\S+)$']
    ['^LoggingEndDate=+' '(\S+)$'] % this variable doesn`t exist in RBRConcerto files
    ['^LoggingEndTime=+' '(\S+\s?\S+)$']
    ['^LoggingSamplingPeriod=+' '(\d+)Hz'] % this variable doesn`t exist in RBRConcerto files
    ['^LoggingSamplingPeriod=+' '(\d\d:\d\d:\d\d)']
    ['^NumberOfChannels=+' '(\d+)']
    ['^CorrectionToConductivity=+' '(\d+)'] % this variable doesn`t exist in RBRConcerto files
    ['^NumberOfSamples=+' '(\d+)']
    ['^HostVersion=+' '(\S+\s?\S+\s?\S+\s?\S+\s?\S+\s?)$'] % this is used to determine date formating of data 5 characters in 1.13.10 legacy export txt
  };
  
  startDate = '';
  startTime = '';
  endDate = '';
  endTime = '';

  for k = 1:length(lines)
    
    % try exprs until we get a match
    for m = 1:length(exprs)
    
      % check for the line containing start sample time
      tkns = regexp(lines{k}, exprs{m}, 'tokens');
      
      if isempty(tkns), continue; end
      
      header.make     = 'RBR';
      
      switch m
          % instrument information
          case 1
              header.model    = tkns{1}{1};
          case 2
              header.firmware = tkns{1}{1};
          case 3
              header.serial   = tkns{1}{1};
              
          % start of sampling
          case 4
              startDate    = tkns{1}{1}; % this variable doesn`t exist in RBRConcerto files
          case 5
              startTime    = tkns{1}{1};
              
          % end of sampling
          case 6
              endDate    = tkns{1}{1}; % this variable doesn`t exist in RBRConcerto files
          case 7
              endTime    = tkns{1}{1};
              
          % sample interval
          case 8 % this variable doesn`t exist in RBRConcerto files
              header.interval = 1/str2double(tkns{1}{1});
              
          % other sample interval
          case 9
              [~, ~, ~, H, MN, S] = datevec(datenum(tkns{1}{1}, 'HH:MM:SS'));
              header.interval = H*3600 + MN*60 + S;
              
          % number of channels
          case 10
              header.channels = str2double(tkns{1}{1});
          
          % correction to conductivity
          case 11 % this variable doesn`t exist in RBRConcerto files
              header.correction  = tkns{1}{1};
              
          % number of samples
          case 12
              header.samples  = str2double(tkns{1}{1});
              
          % Ruskin Version - used to select date format
          case 13
              header.hostversion  = tkns{1}{1};
      end
    end
  end
  % section meant to remedy the changing versions of Ruskin and the
  % nomenclature surrounding start and end dates
  ruskinVer = contains(header.hostversion, {'1.13.7', '1.13.10', '1.13.13', '2.3.0', '2.9.2'}); % to accomodate the latest 2019 version of ruskin
  
  if isempty(ruskinVer)
      ruskinVer = contains(header.hostversion, '1.13.7'); end
      
   if ~isempty(startDate) && ~isempty(startTime) % ruskin v1.5
      if length(startDate) == 8 
          header.start    = datenum([startDate ' ' startTime],  'yy/mm/dd HH:MM:SS.FFF');
      else
          header.start    = datenum([startDate ' ' startTime],  'yyyy/mm/dd HH:MM:SS.FFF');
      end
  end    
  
  if isempty(startDate) && ~isempty(startTime) % ruskin v1.7+ % doesn`t work for V1.13.7 but it does for v1.13.10
     
      if length(startTime) == 11 && ruskinVer == 0
          header.start    = datenum([startDate ' ' startTime],  'dd-mmm-yyyy HH:MM:SS.FFF');
      else
          header.start    = datenum([startDate ' ' startTime],  'yyyy-mm-dd HH:MM:SS.FFF');  
      end
  end
  
  if ~isempty(endDate) && ~isempty(endTime) % ruskin v1.5
      if length(endDate) == 8
          header.end      = datenum([endDate ' ' endTime],    'yy/mm/dd HH:MM:SS.FFF');
      else
          header.end      = datenum([endDate ' ' endTime],    'yyyy/mmm/dd HH:MM:SS.FFF');
      end
  end   
  
  if isempty(endDate) && ~isempty(endTime) % ruskin v1.7+
      if length(endTime) == 11 && ruskinVer == 0
          header.end    = datenum([endDate ' ' endTime],  'dd-mmm-yyyy HH:MM:SS.FFF');
      else
          header.end    = datenum([endDate ' ' endTime],  'yyyy-mm-dd HH:MM:SS.FFF');
      end
  end
  
end

function data = readData(fid, header)
%READDATA Reads the sample data from the file.

  data = struct;
 
  %firmwareNum to be used to select for date- time formatting unique to
  %various versions of Ruskin and firmware
  firmwareNum = str2double(header.firmware);
  
  % section meant to remedy the changing versions of Ruskin and the
  % nomenclature surrounding start and end dates
  ruskinVer = contains(header.hostversion, {'1.13.10', '1.13.13', '2.3.0', '2.9.2'}); % updated to reflect the 2019 ruskin update
  
  if isempty(ruskinVer) % was ~ , shawn sept 18, 2018
      ruskinVer = strfind(header.hostversion, '1.13.7'); end
  
  % get the column names
  header.variables = strrep(header.variables, ' & ', '|');
  header.variables = strrep(header.variables, '  ', '|');
  while ~strcmpi(header.variables, strrep(header.variables, '||', '|'))
      header.variables = strrep(header.variables, '||', '|');
  end
  cols = textscan(header.variables, '%s', 'Delimiter', '|');
  cols = cols{1};
  
  % rename variables with '-', ' ', '&', '(', ')' as Matlab doesn't allow 
  % them within a structure name
  cols = strrep(cols, '-', '');
  cols = strrep(cols, ' ', '');
  cols = strrep(cols, '(', '');
  cols = strrep(cols, ')', '');
  cols = strrep(cols, '&', '');
  
  % first 2 columns are date and time
  fmt  = '%s %s';
  
  % figure out number of columns from the number of channels
  fmt = [fmt repmat(' %f', [1, length(cols)-1])];
  
  % read in the sample data
  samples = textscan(fid, fmt, 'treatAsEmpty', {'null'});
  
  for k = 1:length(cols)
      % check that all columns have the same length. If not correct it.
      if k>1
          while length(samples{k}) < lenData
              samples{k}(end+1) = NaN;
          end
      else
          lenData = length(samples{k});
      end
      
      % save sample data into the data struct, 
      % using  column names as struct field names
      data.(cols{k}) = samples{k}; 
  end
  
   if length(data.Date{1}) == 8 
      data.time = datenum(data.Date, 'yy/mm/dd') + datenum(data.Time, 'HH:MM:SS.FFF') - datenum('00:00:00', 'HH:MM:SS');
  
  elseif isempty(strfind(data.Date{1}, '-'))
          data.time = datenum(data.Date, 'yyyy/mm/dd') + datenum(data.Time, 'HH:MM:SS.FFF') - datenum('00:00:00', 'HH:MM:SS');
          
  % select date-time format by ruskin and firmware version
  elseif firmwareNum == 10.550
          data.time = datenum(data.Date, 'dd-mmm-yyyy') + datenum(data.Time, 'HH:MM:SS.FFF') - datenum('00:00:00', 'HH:MM:SS'); %ruskin v1.12.6
  elseif firmwareNum < 12.1 & firmwareNum > 11 
		  data.time = datenum(data.Date, 'yyyy-mm-dd') + datenum(data.Time, 'HH:MM:SS.FFF') - datenum('00:00:00', 'HH:MM:SS'); % ruskin v1.8.10
  elseif firmwareNum < 7 & ruskinVer > 0
		  data.time = datenum(data.Date, 'yyyy-mm-dd') + datenum(data.Time, 'HH:MM:SS.FFF') - datenum('00:00:00', 'HH:MM:SS'); % ruskin v1.13.7 and 1.13.10
  else 
		  data.time = datenum(data.Date, 'dd-mmm-yyyy') + datenum(data.Time, 'HH:MM:SS.FFF') - datenum('00:00:00', 'HH:MM:SS'); %ruskin v1.11.1
 end
  data = rmfield(data, 'Date');
  data = rmfield(data, 'Time');
end

function data = readDataCsv(filename)
%READDATA Reads the sample data from the CSV file data type, code from ipsRead, possible use here

  data = struct;
    
  fid = fopen(filename, 'rt');
  params = textscan(fid,'%s',3,'Delimiter',',');
  values = textscan(fid,repmat('%f',[1,3]),'Delimiter',',');
  fclose(fid);

% Example of DraftFile CSV type - exported from IPS5Extract 5.0  
%Date[yyyy/mm/dd HH:MM:SS.FFF],Draft(m), DraftError(m)
% 2020/9/20 - 15:52.18.01 , -0.15893, 0.918689

%Y=values{1};
DTime=values{1};
%M=values{2};
%D=values{3};
%H=values{4};
%MN=values{5};
%S=values{6};
Draft = values{2}; % draft values
Draft_err = values{3}; % Draft Error values

% adding DateTime to values and params
%params{1,1}{9,1} = 'DATETIME';
%values{9} = datenum(Y,M,D,H,MN,S); % DATETIME
%DTime = values{9};

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
clear params values Sort_order inRangeA DTime Draft Draft_err c

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
csvfiles=filenamesCSV(endsWith(filenamesCSV,'_data.txt'));
csvfolders=filefolders(endsWith(filenamesCSV,'_data.txt')); 
%
%Make a cell array of strings containing the full file locations of the
%files.
files=fullfile(csvfolders(:),csvfiles(:));
%
%Import the csv tables into out table
% setting the variables to extract
%
for i = 1:length(files)
    %opts = detectImportOptions(files{i},'PreserveVariableNames', true,'VariableNamesLine',1);
    %out{i} = readtable(files{i},opts);
    % above code causes issues in 2017b, so changed the code to below -
    % shawn aug 22, 2020
    opts = detectImportOptions(files{i},'ReadVariableNames',true,'VariableNamesLine',1,'ExtraColumnsRule','ignore', 'Delimiter', ',','ConsecutiveDelimitersRule','join');
    opts.SelectedVariableNames = {'Date_yyyy_mm_ddHH_MM_SS_FFF_','TiltX_deg_','TiltY_deg_','ParosTemperature_C_','Pressure_dbar_','Temperature_C_','Battery_V_'};
    out{i} = readtable(files{i},opts);
    %
    % maybe import col 1 applying a datenum conversion from str to ISO date
end
%
% merge the outputs into a usable table
values = vertcat(out{:});
%  
% Sorting the dateTime of the table in ascending order without duplicates
% new date table is created, though dates are strings 
%
isodates(:,1) = datenum(values{:,1},'yyyy/mm/dd-HH:MM:SS');
%
% new array values with mean of sorted dateTime
for j = 2:7
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