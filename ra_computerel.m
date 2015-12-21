function RELout = ra_computerel(varargin)
%Prepare and execute cmdstan code for dependability analyses
%
%Required Input:
% data - data table outputted from the ra_loadfile script (see ra_loadfile
%  for more information about table format)
%
%Outputs:
%
%
%
%
%
%
%History
% by Peter Clayson (12/15/15)
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
        'See help ra_computerel for more information on optional inputs'));
    end
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('data',varargin),1);
    if ~isempty(ind)
        datatable = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
   
elseif ~isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Optional inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help ra_computerel for more information on optional inputs'));
    
end %if ~isempty(varargin)

%ensure the necessary columns are present in the table
colnames = datatable.Properties.VariableNames;

if ~sum(strcmpi(colnames,'id')) 
    if ~exist('headererror','var')
        headerror{1} = 'Subject ID';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Subject ID';
    end
end

if ~sum(strcmpi(colnames,'meas')) 
    if ~exist('headererror','var')
        headerror{1} = 'Measurement';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'Measurement';
    end
end

if ~sum(strcmpi(colnames,'id2')) 
    if ~exist('headererror','var')
        headerror{1} = 'id2';
    elseif ~exist('headererror','var')
        headerror{end+1} = 'id2';
    end
end

%error catch in case headers for the columns needed are not specified
if exist('headerror','var')
    error('varargin:colheaders',... %Error code and associated error
    strcat('WARNING: Column headers not properly defined in data \n\n',... 
    'Please properly specify the headers for\n',...
    char(strjoin(headerror,', ')),'\n',...
    'See help ra_loadfile for data table format \n'));
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

%determine how cmstan will be set up
%analysis variable will indicate whether group or events need to be
%considered when sending code to cmdstan
%analysis:
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible gruops and event types to consider

if ~exist('ngroup','var') && ~exist('nevent','var')
    analysis = 1;
elseif exist('ngroup','var') && ~exist('nevent','var')
    analysis = 2;
elseif ~exist('ngroup','var') && exist('nevent','var')
    analysis = 3;
elseif exist('ngroup','var') && exist('nevent','var')
    analysis = 4;
end

%settings for cmdstan
niter = 100;
nchains = 3;

%define a dummy variable to store the character for a single quote to be
%used in the code sent to cmdstan
q = char(39);

%create a structure array to store information
REL = struct;
REL.niter = niter;
REL.nchains = nchains;
REL.data = datatable;

switch analysis
    case 1 %no groups or event types to consider
        
        REL.groups = 'none';
        REL.events = 'none';
        
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

        data = struct(...
            'NOBS',length(datatable.id), ... %number of observations
            'NSUB', length(unique(datatable.id)),... %number of participants
            'id', datatable.id2,... %id variable
            'meas', datatable.meas); %measurement variable

        fit = stan('model_code', stan_in, 'model_name', 'test1',...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', true, 'file_overwrite', true);

        fit.block();
        print(fit);    

        REL.out.mu = fit.extract('pars','mu').mu;
        REL.out.sig_u = fit.extract('pars','sig_u').sig_u;
        REL.out.sig_e = fit.extract('pars','sig_e').sig_e; 

    case 2 %possible multiple groups 

        REL.groups = groupnames;
        REL.events = 'none';

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
        
        REL.stan_in = stan_in;
        
        %create structure array for the data to be sent to cmdstan
        data = struct;
        
        for i=1:ngroup
            fieldname = sprintf('NG%d',i);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end 
        
        data = struct('NDEP',length(erndep.id), ... %number of observations
    'NANX',length(ernanx.id), ... %number of observations
    'NCNT',length(erncnt.id), ... %number of observations
    'JDEP', length(unique(erndep.id)),... %number of participants
    'JANX', length(unique(ernanx.id)),... %number of participants
    'JCNT', length(unique(erncnt.id)),... %number of participants
    'id_dep', erndep.id2,... %id variable
    'id_anx', ernanx.id2,... %id variable
    'id_cnt', erncnt.id2,... %id variable
    'roi_dep', erndep.ern_roi,...
    'roi_anx', ernanx.ern_roi,...
    'roi_cnt', erncnt.ern_roi); %measurement variable

      
        fit = stan('model_code', stan_in, 'model_name', 'test1',...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', true, 'file_overwrite', true);

        fit.block();
        print(fit);    

    case 3 %possible event types but no groups to consider



    case 4 %possible gruops and event types to consider





end





end

