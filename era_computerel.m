function RELout = era_computerel(varargin)
%Prepare and execute cmdstan code for dependability analyses
%
%era_computerel('data',era_datatable,'chains',3,'iter',1000)
%
%Lasted Updated 4/27/16
%
%Required Inputs:
% data - data table outputted from the era_loadfile script (see era_loadfile
%  for more information about table format)
%  Note: era_loadfile sets up the datatable in a specific format with
%  specific header names that are used in this script
% chains - number of chains to run in stan 
% iter - number of iterations to run in stan 
%
%Optional Inputs:
% verbose - 1: Do not print iterations, 2: Print iterations
%
%Outputs:
% RELout - structure array with the following fields.
%  filename: filename of the processed dataset
%  niter: number of interations run in Stan
%  events: names of the events (if applicable)
%  groups: names of the groups (if applicable)
%  data: compiled tables that were used to create data structures to pass
%   to stan
%  out: data output from stan and labels corresponding to order of data
%   mu: person mean
%   sig_u: person variance
%   sig_e: error variance
%   labels: labels for the data. If either events or groups were processed
%    then this the labeles for events or groups. If both events and groups
%    were processed then this labels will be 'event_group' for each
%    possible combination (e.g., 'Error_;_Controls' or 'Hits_;_Patients'  
%    event and group will be separated by _;_ to avoid using a common 
%    delimeter that could be used by a researcher to label an event)
%  conv: nest for convergence data
%   data: convergence statistics with three columns (parameter, n_eff,
%    r_hat)
%  stan_in: compiled model syntax that was passed to Stan

%For more information about how variance components were estimated please
% see
%
% Baldwin, S. A., Larson, M. J., & Clayson, P. E. (2015). The dependability 
% of electrophysiological measurements of performance monitoring in a 
% clinical sample: A generalizability and decision analysis of the ERN and 
% Pe. Psychophysiology, 52(6), 790-800. http://doi.org/10.1111/psyp.12401
%
%It would be terribly remiss of me not thank Dr. Scott Baldwin for
%conceptualizing and developing the formulas that are implemented in this
%toolbox. Dr. Baldwin also wrote the original Stan syntax in R and 
%graciously provided me with all of his code. This Matlab code is based
%off of his R code. 

% Copyright (C) 2016 Peter E. Clayson
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program (gpl.txt). If not, see 
%     <http://www.gnu.org/licenses/>.
%

%History
% by Peter Clayson (4/18/16)
% peter.clayson@gmail.com

%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_computerel for more information on optional inputs'));
    end
    
    %check if the dataset is present
    ind = find(strcmp('data',varargin),1);
    if ~isempty(ind)
        datatable = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
    
    %check whether chains is specified
    ind = find(strcmp('chains',varargin),1);
    if ~isempty(ind)
        nchains = varargin{ind+1}; 
    else 
        error('varargin:nchains',... %Error code and associated error
        strcat('WARNING: Number of chains not specified \n\n',... 
        'Please input the number of chains for stan\n',...
        'See help era_computerel for more information on inputs'));
    end
    
    %check whether iter is specified
    ind = find(strcmp('iter',varargin),1);
    if ~isempty(ind)
        niter = varargin{ind+1}; 
    else 
        error('varargin:niter',... %Error code and associated error
        strcat('WARNING: Number of iterations not specified \n\n',... 
        'Please input the number of iterations for stan\n',...
        'See help era_computerel for more information on inputs'));
    end
    
    %check whether verbose is specified
    ind = find(strcmp('verbose',varargin),1);
    if ~isempty(ind)
        verbose = varargin{ind+1}; 
    else 
        verbose = 1;
    end
    
elseif isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help era_computerel for more information on inputs'));
    
end %if ~isempty(varargin)

eraver = '0.3.2';

%store raw data to pass to era_checkconvergence in the event that the user
%chooses to run with more iterations
dataraw = datatable;

%ensure the necessary columns are present in the table (at least the
%headers for id and meas)
colnames = datatable.Properties.VariableNames;

%check for id column
if ~sum(strcmpi(colnames,'id')) 
    if ~exist('headererror','var')
        headerror{1} = 'Subject ID';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Subject ID';
    end
end

%check for measurement column
if ~sum(strcmpi(colnames,'meas')) 
    if ~exist('headererror','var')
        headerror{1} = 'Measurement';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Measurement';
    end
end

%error catch in case headers for the columns needed are not specified
if exist('headerror','var')
    error('varargin:colheaders',... %Error code and associated error
    strcat('WARNING: Column headers not properly defined in data \n\n',... 
    'Please properly specify the headers for\n',...
    char(strjoin(headerror,', ')),'\n',...
    'See help era_loadfile for data table format \n'));
end

%change verbosity to match true or false
if verbose == 1
    verbosity = false;
elseif verbose == 2
    verbosity = true;
end

%check whether groups or events are in the table
if sum(strcmpi(colnames,'group'))
    groupnames = unique(datatable.group(:));
    ngroup = length(groupnames);
end

if sum(strcmpi(colnames,'event'))
    eventnames = unique(datatable.event(:));
    nevent = length(eventnames);
end

%create a structure array to store information that will later be outputted
REL = struct;
REL.filename = datatable.Properties.Description;
REL.niter = niter;
REL.nchains = nchains;

%determine how cmdstan will be set up
%analysis variable will indicate whether group or events need to be
%considered when sending code to cmdstan
%analysis:
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

if ~exist('ngroup','var') && ~exist('nevent','var')
    analysis = 1;
elseif exist('ngroup','var') && ~exist('nevent','var')
    analysis = 2;
elseif ~exist('ngroup','var') && exist('nevent','var')
    analysis = 3;
elseif exist('ngroup','var') && exist('nevent','var')
    analysis = 4;
end

switch analysis
    case 1 %no groups or event types to consider
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,'id');
        id2 = zeros(0,length(datatable.id));

        for i = 1:length(datatable.id)
            if i == 1
                id2(1) = 1;
                count = 1;
            elseif i > 1 && strcmp(char(datatable.id(i)), char(datatable.id(i-1)))
                id2(i) = count;
            elseif i > 1 && ~strcmp(datatable.id(i), datatable.id(i-1))
                count = count+1;
                id2(i) = count;
            end
        end

        datatable.id2 = id2(:);
        
        %store the datatable used in REL
        REL.data = datatable;
        
        %groups and events were not considered
        REL.groups = 'none';
        REL.events = 'none';
        
        %create cmdstand syntax
        stan_in = {
          'data {' 
          '  int<lower=0> NOBS; //number of obs'
          '  int<lower=0> NSUB; //number of subj'
          '  int<lower=0, upper=NSUB> id[NOBS]; //subject id;'
          '  real meas[NOBS];'
          '}'
          'parameters {'
          '  real mu;'
          '  real<lower=0> sig_u;'
          '  real<lower=0> sig_e;'
          '  vector[NSUB] u_raw;'
          '}'
          'transformed parameters {'
          '  vector[NSUB] u;'
          '  u <- mu + sig_u*u_raw;'
          '}'
          'model {'
          '  u_raw ~ normal(0,1);'
          '  for (i in 1:NOBS) {'
          '    meas[i] ~ normal(u[id[i]], sig_e);'
          '  }'
          '  '
          '  mu ~ normal(0,100);'
          '  sig_u ~ cauchy(0,40);'
          '  sig_e ~ cauchy(0,40);'
          '}'
        };
        
        %create data structure to pass to cmdstan
        data = struct(...
            'NOBS',length(datatable.id), ... %number of observations
            'NSUB', length(unique(datatable.id)),... %number of participants
            'id', datatable.id2,... %id variable
            'meas', datatable.meas); %measurement variable
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        modelname = strcat('cmdstan',char(date));
        
        %fit model in cmdstan
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', verbosity, 'file_overwrite', true);
        
        %don't let the user continue to use the Matlab command window
        fit.block();
        
        %extract and store cmdstan output
        REL.out.mu = fit.extract('pars','mu').mu;
        REL.out.sig_u = fit.extract('pars','sig_u').sig_u;
        REL.out.sig_e = fit.extract('pars','sig_e').sig_e; 
        
        %label is simply measure (no events or groups to deal with)
        REL.out.labels = 'measure';
        
        REL.out.conv.data = era_storeconv(fit);
        
        
    case 2 %possible multiple groups 

        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'group','id'});
        id2 = zeros(0,length(datatable.id));

        for i = 1:length(datatable.id)
            if i == 1
                id2(1) = 1;
                count = 1;
            elseif i > 1 &&...
                    strcmp(char(datatable.id(i)), char(datatable.id(i-1)))...
                    &&...
                    strcmp(char(datatable.group(i)), char(datatable.group(i-1)))
                id2(i) = count;
            elseif i > 1 &&...
                    ~strcmp(char(datatable.id(i)), char(datatable.id(i-1)))...
                    &&...
                    strcmp(char(datatable.group(i)), char(datatable.group(i-1)))
                count = count+1;
                id2(i) = count;
            elseif i > 1 &&...
                    ~strcmp(char(datatable.id(i)), char(datatable.id(i-1)))...
                    &&...
                    ~strcmp(char(datatable.group(i)), char(datatable.group(i-1)))
                count = 1;
                id2(i) = 1;
            end
        end

        datatable.id2 = id2(:);
        
        %create labels for the gorups
        REL.data = datatable;
        groupnames = unique(datatable.group(:));
        ngroup = length(groupnames);
        
        %all names are forced to be strings for ease of use later on
        if isnumeric(groupnames)
            groupnames = num2str(groupnames);
        end
        
        %store group names. separate events are not a consideration
        REL.groups = groupnames;
        REL.events = 'none';
        
        
        %create model syntax for cmdstan
        stan_in{1,1} = 'data {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> NG%d; //number of obs in g%d',i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> JG%d; //number of subj in g%d',i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0, upper=JG%d> id_G%d[NG%d]; //subject id g%d;',...
                i,i,i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real meas_G%d[NG%d];',i,i);
        end    

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'parameters {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real mu_G%d;',i);
        end 

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_u_G%d;',i);
        end 

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_e_G%d;',i);
        end 

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  vector[JG%d] u_raw_G%d;',i,i);
        end 

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'transformed parameters {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  vector[JG%d] u_G%d;',i,i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  u_G%d <- mu_G%d + sig_u_G%d*u_raw_G%d;',...
                i,i,i,i);
        end

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'model {';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  u_raw_G%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  for (i in 1:NG%d) {;',i);
            stan_in{end+1,1} = ...
                sprintf('    meas_G%d[i] ~ normal(u_G%d[id_G%d[i]], sig_e_G%d);',...
                i,i,i,i);
            stan_in{end+1,1} = '  }';
        end

        stan_in{end+1,1} = '  ';

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  mu_G%d ~ normal(0,100);',i);
        end

        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  sig_u_G%d ~ cauchy(0,40);',i);
        end
        
        for i=1:ngroup
            stan_in{end+1,1} = ...
                sprintf('  sig_e_G%d ~ cauchy(0,40);',i);
        end
        
        stan_in{end+1,1} = '}';
        
        
        %store cmdstan model syntax
        REL.stan_in = stan_in;
        
        %create structure array for the data to be sent to cmdstan
        data = struct;
        
        for i=1:ngroup
            fieldname = sprintf('NG%d',i);
            fieldvalue = height(datatable(ismember(datatable.group,...
                groupnames(i)),1)); 
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:ngroup
            fieldname = sprintf('JG%d',i);
            dummy = datatable(ismember(datatable.group,...
                groupnames(i)),:);
            fieldvalue = length(unique(dummy.id));
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:ngroup
            fieldname = sprintf('id_G%d',i);
            fieldvalue = datatable{ismember(datatable.group,...
                groupnames(i)),{'id2'}};
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:ngroup
            fieldname = sprintf('meas_G%d',i);
            fieldvalue = datatable{ismember(datatable.group,...
                groupnames(i)),{'meas'}};
            data.(fieldname) = fieldvalue;
        end
        
        %execute model code 
        
        modelname = strcat('cmdstan',char(date));
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        %run cmdstan
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', verbosity, 'file_overwrite', true);
        
        %block user from using Matlab command window
        fit.block();
        
        %create fields for storing parsed cmdstan outputs
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};  
        
        %store cmdstand outputs
        for i=1:ngroup
            measname = sprintf('mu_G%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.mu(:,end+1) = measvalue;           
        end
        
        for i=1:ngroup
            measname = sprintf('sig_u_G%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_u(:,end+1) = measvalue;           
        end
        
        for i=1:ngroup
            measname = sprintf('sig_e_G%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_e(:,end+1) = measvalue;           
        end
        
        %store labels
        for i=1:ngroup
            REL.out.labels(:,end+1) = groupnames(i);           
        end
        
        REL.out.conv.data = era_storeconv(fit);
        
    case 3 %possible event types but no groups to consider

        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'id','event'}); 
        
        eventnames = unique(datatable.event(:));
        nevent = length(eventnames);
        
        %ensure all event names are strings (makes life easier later)
        if isnumeric(eventnames)
            eventnames = num2str(eventnames);
        end
        
        %store the event names
        REL.events = eventnames;
        
        %no groups to consider
        REL.groups = 'none';
        
        %create data structure for storing measurements for each event
        eventarray = struct;
        eventarray.data = [];
        
        for i = 1:nevent
            
            dummytable = datatable(ismember(datatable.event,...
                char(eventnames(i))),:);
            eventarray.data{end+1} = dummytable;
            
        end
        
        %creating that id2 variable to make cmdstan happy
        for j = 1:nevent
            
            id2 = zeros(0,height(eventarray.data{j}));
            
            for i = 1:height(eventarray.data{j})
                
                if i == 1
                    id2(1) = 1;
                    count = 1;
                elseif i > 1 && strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1)))
                    id2(i) = count;
                elseif i > 1 && ~strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1)))
                    count = count+1;
                    id2(i) = count;
                end
            end
            
            eventarray.data{j}.id2 = id2';
            
        end
        
        %store the data
        REL.data = eventarray.data;

        stan_in{1,1} = 'data {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> NE%d; //number of obs in E%d',i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0> JE%d; //number of subj in E%d',i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  int<lower=0, upper=JE%d> id_E%d[NE%d]; //subject id E%d;',...
                i,i,i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real meas_E%d[NE%d];',i,i);
        end    

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'parameters {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real mu_E%d;',i);
        end 

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_u_E%d;',i);
        end 

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  real<lower=0> sig_e_E%d;',i);
        end 

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  vector[JE%d] u_raw_E%d;',i,i);
        end 

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'transformed parameters {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  vector[JE%d] u_E%d;',i,i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  u_E%d <- mu_E%d + sig_u_E%d*u_raw_E%d;',...
                i,i,i,i);
        end

        stan_in{end+1,1} = '}';
        stan_in{end+1,1} = 'model {';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  u_raw_E%d ~ normal(0,1);',i);
            stan_in{end+1,1} = ...
                sprintf('  for (i in 1:NE%d) {;',i);
            stan_in{end+1,1} = ...
                sprintf('    meas_E%d[i] ~ normal(u_E%d[id_E%d[i]], sig_e_E%d);',...
                i,i,i,i);
            stan_in{end+1,1} = '  }';
        end

        stan_in{end+1,1} = '  ';

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  mu_E%d ~ normal(0,100);',i);
        end

        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  sig_u_E%d ~ cauchy(0,40);',i);
        end
        
        for i=1:nevent
            stan_in{end+1,1} = ...
                sprintf('  sig_e_E%d ~ cauchy(0,40);',i);
        end
        
        stan_in{end+1,1} = '}';
        
        REL.stan_in = stan_in;
        
        %create structure array for the data to be sent to cmdstan
        data = struct;

        for i=1:nevent %number of observations
            fieldname = sprintf('NE%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = height(dummytable); 
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:nevent %number of participants
            fieldname = sprintf('JE%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = length(unique(dummytable.id));
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:nevent %id variable
            fieldname = sprintf('id_E%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = dummytable.id2;
            data.(fieldname) = fieldvalue;
        end
        
        for i=1:nevent %measurement variable
            fieldname = sprintf('meas_E%d',i);
            dummytable = eventarray.data{i};
            fieldvalue = dummytable.meas;
            data.(fieldname) = fieldvalue;
        end
        
        %execute model code 
        
        modelname = strcat('cmdstan',char(date));
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        %run in cmdstan
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', verbosity, 'file_overwrite', true);
        
        %block user from using Matlab command window 
        fit.block();
    
        
        %create fields to store cmdstan output
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};
        REL.out.conv.data = {};
        
        %parse outputs
        for i=1:nevent
            measname = sprintf('mu_E%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.mu(:,end+1) = measvalue;           
        end
        
        for i=1:nevent
            measname = sprintf('sig_u_E%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_u(:,end+1) = measvalue;           
        end
        
        for i=1:nevent
            measname = sprintf('sig_e_E%d',i);
            measvalue = fit.extract('pars',measname).(measname);
            REL.out.sig_e(:,end+1) = measvalue;           
        end
       
        %labels are the event names
        for i=1:nevent
            REL.out.labels(:,end+1) = eventnames(i);           
        end
       
        REL.out.conv.data = era_storeconv(fit);
        
    case 4 %possible groups and event types to consider

        %to decrease memory load, cmdstan will be run separately for each
        %event. Groups will be processed in the same cmdstan run.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'group','id','event'});
        eventnames = unique(datatable.event(:));
        if isnumeric(eventnames)
            eventnames = num2str(eventnames);
        end
        
        %store event names
        nevent = length(eventnames);
        REL.events = eventnames;
       
        groupnames = unique(datatable.group(:));
        if isnumeric(groupnames)
            groupnames = num2str(groupnames);
        end
        
        %store group names
        ngroup = length(groupnames);
        REL.groups = groupnames;
        
        %create structure for storing data
        eventarray = struct;
        eventarray.data = [];
        
        %store the data for each event separately
        for i = 1:nevent
            
            dummytable = datatable(ismember(datatable.event,...
                char(eventnames(i))),:);
            eventarray.data{end+1} = dummytable;
            
        end
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement. This is
        %done separately for each event since the event data will be passed
        %separately to cmdstan
        
        for j = 1:nevent
            
            id2 = zeros(0,height(eventarray.data{j}));
            
            for i = 1:height(eventarray.data{j})
                if i == 1
                    id2(1) = 1;
                    count = 1;
                elseif i > 1 &&...
                        strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1))) &&...
                        strcmp(char(eventarray.data{j}.group(i)),...
                        char(eventarray.data{j}.group(i-1)))
                    id2(i) = count;
                elseif i > 1 &&...
                        ~strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1))) &&...
                        strcmp(char(eventarray.data{j}.group(i)),...
                        char(eventarray.data{j}.group(i-1)))
                    count = count+1;
                    id2(i) = count;
                elseif i > 1 &&...
                        ~strcmp(char(eventarray.data{j}.id(i)),...
                        char(eventarray.data{j}.id(i-1))) &&...
                        ~strcmp(char(eventarray.data{j}.group(i)),...
                        char(eventarray.data{j}.group(i-1)))
                    count = 1;
                    id2(i) = 1;
                end
            end
            
            eventarray.data{j}.id2 = id2';
            
        end
        
        %put the compiled data into REL structure
        REL.data = eventarray.data;
        
        %create the fields for storing the parsed cmdstan outputs 
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};
        REL.out.conv.data = {};
        REL.stan_in = {};
        
        %separate syntax will need to be generated for each event because
        %the number of samples will be different
        
        %create syntax for cmdstan
        for j=1:nevent
            
            clear stan_in
            stan_in{1,1} = 'data {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  int<lower=0> NG%d; //number of obs in g%d',i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  int<lower=0> JG%d; //number of subj in g%d',i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  int<lower=0, upper=JG%d> id_G%d[NG%d]; //subject id g%d;',...
                    i,i,i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real meas_G%d[NG%d];',i,i);
            end    

            stan_in{end+1,1} = '}';
            stan_in{end+1,1} = 'parameters {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real mu_G%d;',i);
            end 

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real<lower=0> sig_u_G%d;',i);
            end 

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  real<lower=0> sig_e_G%d;',i);
            end 

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  vector[JG%d] u_raw_G%d;',i,i);
            end 

            stan_in{end+1,1} = '}';
            stan_in{end+1,1} = 'transformed parameters {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  vector[JG%d] u_G%d;',i,i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  u_G%d <- mu_G%d + sig_u_G%d*u_raw_G%d;',...
                    i,i,i,i);
            end

            stan_in{end+1,1} = '}';
            stan_in{end+1,1} = 'model {';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  u_raw_G%d ~ normal(0,1);',i);
                stan_in{end+1,1} = ...
                    sprintf('  for (i in 1:NG%d) {;',i);
                stan_in{end+1,1} = ...
                    sprintf('    meas_G%d[i] ~ normal(u_G%d[id_G%d[i]], sig_e_G%d);',...
                    i,i,i,i);
                stan_in{end+1,1} = '  }';
            end

            stan_in{end+1,1} = '  ';

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  mu_G%d ~ normal(0,100);',i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  sig_u_G%d ~ cauchy(0,40);',i);
            end

            for i=1:ngroup
                stan_in{end+1,1} = ...
                    sprintf('  sig_e_G%d ~ cauchy(0,40);',i);
            end

            stan_in{end+1,1} = '}';
            
            %store the cmdstan syntax
            REL.stan_in{end+1} = stan_in;

            %create structure for the data to be sent to cmdstan
            data = struct;

            %prep data structure for cmdstan
            for i=1:ngroup
                fieldname = sprintf('NG%d',i);
                fieldvalue = length(eventarray.data{j}.group(ismember(eventarray.data{j}.group,...
                    groupnames(i)),1)); 
                data.(fieldname) = fieldvalue;
            end

            for i=1:ngroup
                fieldname = sprintf('JG%d',i);
                dummy = eventarray.data{j}(ismember(eventarray.data{j}.group,...
                    groupnames(i)),:);
                fieldvalue = length(unique(dummy.id));
                data.(fieldname) = fieldvalue;
            end

            for i=1:ngroup
                fieldname = sprintf('id_G%d',i);
                fieldvalue = eventarray.data{j}{ismember(eventarray.data{j}.group,...
                    groupnames(i)),{'id2'}};
                data.(fieldname) = fieldvalue;
            end

            for i=1:ngroup
                fieldname = sprintf('meas_G%d',i);
                fieldvalue = eventarray.data{j}{ismember(eventarray.data{j}.group,...
                    groupnames(i)),{'meas'}};
                data.(fieldname) = fieldvalue;
            end

            %execute model code 

            modelname = strcat('cmdstan',char(date));

            fprintf('\nModel is being run in cmdstan\n');
            fprintf('\nThis may take a while depending on the amount of data\n');
            
            %run cmdstan
            fit = stan('model_code', stan_in, 'model_name', modelname,...
                'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
                niter/10, 'verbose', verbosity, 'file_overwrite', true);
            
            %block user from using Matlab command window
            fit.block();
            
            %parse cmdstan outputs
            for i=1:ngroup
                measname = sprintf('mu_G%d',i);
                measvalue = fit.extract('pars',measname).(measname);
                REL.out.mu(:,end+1) = measvalue;           
            end

            for i=1:ngroup
                measname = sprintf('sig_u_G%d',i);
                measvalue = fit.extract('pars',measname).(measname);
                REL.out.sig_u(:,end+1) = measvalue;           
            end

            for i=1:ngroup
                measname = sprintf('sig_e_G%d',i);
                measvalue = fit.extract('pars',measname).(measname);
                REL.out.sig_e(:,end+1) = measvalue;           
            end

            %store labels for each column by separating the event name and
            %group name with an underscore
            for i=1:ngroup
                REL.out.labels(:,end+1) = strcat(eventnames(j),...
                    '_;_',groupnames(i));           
            end

            REL.out.conv.data{end+1} = era_storeconv(fit);
            
        end %for j=1:nevent

end %switch analysis

REL.eraver = eraver;
RELout = REL;

end

function convstats = era_storeconv(fit)
%pull out r_hats from Stanfit model

%pull summary using the print function
%there's no other way to get at the r_hat values, so output will also be
%printed in the command window
output = print(fit);

%define cell array that holds the r_hats
convstats = cell(0,3);
convstats{1,1} = 'name';
convstats{1,2} = 'n_eff';
convstats{1,3} = 'r_hat';

%cycle through the important inputs
%depending on whether there are multiple events/groups there may be 
%multiple mu_, sig_u_, and sig_e_
inp2check = {'lp__' 'mu_' 'sig_u_' 'sig_e_'}; 

for j = 1:length(inp2check)
    check = strncmp(output,inp2check{j},length(inp2check{j}));
    ind2check = find(check == 1);
    for i = 1:length(ind2check)
        pullrow = strsplit(output{ind2check(i),:},' ');
        label = pullrow{1};
        neff = str2num(pullrow{end-2});
        rhat = str2num(pullrow{end});
        convstats{end+1,1} = label;
        convstats{end,2} = neff;
        convstats{end,3} = rhat;
    end
end

end

