function RELout = ra_computerel(varargin)
%Prepare and execute cmdstan code for dependability analyses
%
%Required Input:
% data - data table outputted from the ra_loadfile script (see ra_loadfile
%  for more information about table format)
%
%Outputs:
% RELout - structure array with the following fields.
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

%settings for cmdstan
niter = 100;
nchains = 3;

%create a structure array to store information
REL = struct;
REL.filename = datatable.Properties.Description;
REL.niter = niter;
REL.nchains = nchains;


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

        REL.data = datatable;
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
        
        fprintf('\nModel is being run in cmdstan\n');
        fprintf('\nThis may take a while depending on the amount of data\n');
        
        modelname = strcat('cmdstan',char(date));
        
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', false, 'file_overwrite', true);
        
        fit.block();
        
        REL.stanfit = fit;
        
        REL.out.mu = fit.extract('pars','mu').mu;
        REL.out.sig_u = fit.extract('pars','sig_u').sig_u;
        REL.out.sig_e = fit.extract('pars','sig_e').sig_e; 
        REL.out.labels = 'measure';
        
        
        
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
        
        REL.data = datatable;
        groupnames = unique(datatable.group(:));
        ngroup = length(groupnames);
        
        if isnumeric(groupnames)
            groupnames = num2str(groupnames);
        end
        
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
        
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', false, 'file_overwrite', true);

        fit.block();
        
           

        REL.stanfit = fit;
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};  
        
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
       
        for i=1:ngroup
            REL.out.labels(:,end+1) = groupnames(i);           
        end
        
        
        
    case 3 %possible event types but no groups to consider

        
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'id','event'}); 
        
        eventnames = unique(datatable.event(:));
        nevent = length(eventnames);
        
        if isnumeric(eventnames)
            eventnames = num2str(eventnames);
        end
        
        REL.events = eventnames;
        REL.groups = 'none';
        
        eventarray = struct;
        eventarray.data = [];
        
        for i = 1:nevent
            
            dummytable = datatable(ismember(datatable.event,...
                char(eventnames(i))),:);
            eventarray.data{end+1} = dummytable;
            
        end
        
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
        
        
        REL.data = eventarray.data;
        
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};
        
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
        
        fit = stan('model_code', stan_in, 'model_name', modelname,...
            'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
            niter/10, 'verbose', false, 'file_overwrite', true);
        
        fit.block();
           

        REL.stanfit = fit;  
    
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
       
        for i=1:nevent
            REL.out.labels(:,end+1) = eventnames(i);           
        end
       
        
        
    case 4 %possible groups and event types to consider

        %to decrease memory load, cmdstan will be run separately for each
        %event. Groups will be processed in the same cmdstan run.
        
        %cmdstan requires the id variable to be numeric and sequential. 
        %an id2 variable is created to satisfy this requirement.
        
        fprintf('\nPreparing data for analysis...\n');
        
        datatable = sortrows(datatable,{'group','id','event'});
        eventnames = unique(datatable.event(:));
        if isnumeric(eventnames)
            eventnames = num2str(eventnames);
        end
        
        nevent = length(eventnames);
        REL.events = eventnames;
        
        groupnames = unique(datatable.group(:));
        if isnumeric(groupnames)
            groupnames = num2str(groupnames);
        end
        
        ngroup = length(groupnames);
        REL.groups = groupnames;
        
        eventarray = struct;
        eventarray.data = [];
        
        for i = 1:nevent
            
            dummytable = datatable(ismember(datatable.event,...
                char(eventnames(i))),:);
            eventarray.data{end+1} = dummytable;
            
        end
        
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
        
        REL.data = eventarray.data;
        
        REL.out = [];
        REL.out.mu = [];
        REL.out.sig_u = [];
        REL.out.sig_e = [];
        REL.out.labels = {};
        REL.stanfit = {};
        
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

            REL.stan_in = stan_in;

            %create structure array for the data to be sent to cmdstan
            data = struct;

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
            
            fit = stan('model_code', stan_in, 'model_name', modelname,...
                'data', data, 'iter', niter,'chains', nchains, 'refresh',... 
                niter/10, 'verbose', false, 'file_overwrite', true);
            
            fit.block();
                
            
            REL.stanfit{end+1} = fit; 
            
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

            for i=1:ngroup
                REL.out.labels(:,end+1) = strcat(eventnames(j),...
                    '_',groupnames(i));           
            end

        end

end

RELout = REL;

end

