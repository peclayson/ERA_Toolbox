function RELout = ra_relfigures(varargin)
%Creates figures depicting dependability estimates and displays information
% about optimal cutoffs and overall dependability
%
%ra_relfigures('data',REL)
%
%Note: The Statistics and Machine Learning Toolbox is required
%
%Required Inputs:
% data - structure array containing results of dependability analyses using
%  cmdstan. See ra_computerel for more information.
%
%Optional Inputs:
% relcutoff - reliability level to use for cutoff when deciding the
%  minimum number of trials needed to achieve this specified level of
%  reliability
%
%Output:
% One figure is plotted that displays the dependability of measurements as
%  the number of trials increases. If there is separate information for
%  groups, events, or groups and events, the dependability estimates will
%  be displayed accordingly.

%History 
% by Peter Clayson (12/23/15)
% peter.clayson@gmail.com
%


if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help ra_relfigures for more information on optional inputs'));
    end
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('data',varargin),1);
    if ~isempty(ind)
        REL = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
   
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('relcutoff',varargin),1);
    if ~isempty(ind)
        relcutoff = cell2mat(varargin{ind+1}); 
    else 
        relcutoff = .70; %default level is .70
    end
    
elseif ~isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Optional inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help ra_relfigures for more information on optional inputs'));
    
end %if ~isempty(varargin)

data = struct;

if strcmpi(REL.groups,'none')
    ngroups = 1;
    gnames = cellstr(REL.groups);
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

if strcmpi(REL.events,'none')
    nevents = 1;
    enames = cellstr(REL.events);
else
    nevents = length(REL.events);
    enames = REL.events(:);
end


%figure out whether groups or events need to be considered
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

if ngroups == 1 && nevents == 1
    analysis = 1;
elseif ngroups > 1 && nevents == 1
    analysis = 2;
elseif ngroups == 1 && nevents > 1
    analysis = 3;
elseif ngroups > 1 && nevents > 1
    analysis = 4;
end

switch analysis
    case 1 %1 - no groups or event types to consider

        gloc = 1;
        eloc = 1;
        
        data.g(gloc).e(eloc).label = REL.out.labels(gloc);
        data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
        data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,gloc);
        data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,gloc);
        data.g(gloc).e(eloc).elabel = cellstr('none');
        data.g(gloc).glabel = gnames(gloc);
        
    case 2 %2 - possible multiple groups but no event types to consider 
        
        eloc = 1;
        
        for gloc=1:length(REL.out.labels)
            
            data.g(gloc).e(eloc).label = REL.out.labels(gloc);
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,gloc);
            data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,gloc);
            data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,gloc);
            data.g(gloc).e(eloc).elabel = cellstr('none');
            data.g(gloc).glabel = gnames(gloc);
    
        end
        
    case 3 %3 - possible event types but no groups to consider
        
        gloc = 1;
        
        for eloc=1:length(REL.out.labels)

            data.g(gloc).e(eloc).label = REL.out.labels(eloc);
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,eloc);
            data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,eloc);
            data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,eloc);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
    
        end
        
    case 4 %4 - possible groups and event types to consider
        for i=1:length(REL.out.labels)
    
            lblstr = strsplit(REL.out.labels{i},'_');

            eloc = find(ismember(enames,lblstr(1)));
            gloc = find(ismember(gnames,lblstr(2)));

            if isempty(eloc); eloc = 1; end;
            if isempty(gloc); gloc = 1; end;

            data.g(gloc).e(eloc).label = REL.out.labels(i);
            data.g(gloc).e(eloc).mu.raw = REL.out.mu(:,i);
            data.g(gloc).e(eloc).sig_u.raw = REL.out.sig_u(:,i);
            data.g(gloc).e(eloc).sig_e.raw = REL.out.sig_e(:,i);
            data.g(gloc).e(eloc).elabel = enames(eloc);
            data.g(gloc).glabel = gnames(gloc);
    
        end
end
  


%flexible number of trials to display (warn about extrapolating too far!)
ntrials = 50;
x = 1:ntrials;

mrel = zeros(ntrials,0);

if nevents > 2
    xplots = ceil(sqrt(nevents));
    yplots = ceil(sqrt(nevents));
elseif nevents == 2
    xplots = 1;
    yplots = 2;
elseif nevents == 1
    xplots = 1;
    yplots = 1;
end

fig = figure;
fig.Position = [230 100 1070 810];
fsize = 16;

for eloc=1:nevents
    for gloc=1:ngroups
        for trial=1:ntrials
            mrel(trial,gloc) = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial)); 
        end
    end
    subplot(yplots,xplots,eloc);
    plot(x,mrel);
    axis([0 ntrials 0 1]);
    set(gca,'fontsize',16);
    
    if ~strcmpi(enames{eloc},'none')
        title(enames{eloc},'FontSize',20);
    end
    
    ylabel('Dependability','FontSize',fsize);
    xlabel('Number of Observations','FontSize',fsize);
    hline = refline(0,.7);
    set(hline,'Color','b','LineStyle',':');
    if analysis ~= 1 && analysis ~= 3
        leg = legend(gnames{:},'Location','southeast');
        set(leg,'FontSize',fsize);
    end
end

mrel = zeros(0,ntrials);
llrel = zeros(0,ntrials);
ulrel = zeros(0,ntrials);

relsummary.relcutoff = relcutoff;

switch analysis
    case 1 %no groups or event types to consider
        
        eloc = 1;
        gloc = 1;
        
        relsummary.group(gloc).name = gnames{gloc};
        relsummary.group(gloc).event(eloc).name = 'measure';

        for trial=1:ntrials
            mrel(trial) = ...
                mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial)); 
            llrel(trial) = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial),.025);
            ulrel(trial) = quantile(reliab(...
                data.g(gloc).e(eloc).sig_u.raw,...
                data.g(gloc).e(eloc).sig_e.raw,trial),.975);
        end

        %find the number of trials to reach cutoff based on the
        %lower limit of the confidence interval
        trlcutoff = find(llrel >= relcutoff, 1);

        relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
        relsummary.group(gloc).event(eloc).mrel = mrel(trlcutoff);
        relsummary.group(gloc).event(eloc).llrel = llrel(trlcutoff);
        relsummary.group(gloc).event(eloc).ulrel = ulrel(trlcutoff);


        %Only calculate overall dependability on the ids with
        %enough trials

        datatrls = REL.data;

        trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

        ind2include = trltable.GroupCount >= trlcutoff;
        ind2exclude = trltable.GroupCount < trlcutoff;

        relsummary.group(gloc).event(eloc).eventgoodids =...
            trltable.id(ind2include);
        relsummary.group(gloc).event(eloc).eventbadids =...
            trltable.id(ind2exclude);
        
        relsummary.group(gloc).goodids = ...
            relsummary.group(gloc).event(eloc).eventgoodids;
        relsummary.group(gloc).badids = ...
            relsummary.group(gloc).event(eloc).eventbadids;

        datatable = REL.data;

        goodids = table(relsummary.group(gloc).goodids);

        lmedata = innerjoin(datatable, goodids,...
            'LeftKeys', 'id', 'RightKeys', 'Var1',...
            'LeftVariables', {'id' 'meas'});

        trltable = varfun(@length,lmedata,...
            'GroupingVariables',{'id'});

        trlmean = mean(trltable.GroupCount);

        %perform calculations for overall reliability
        lmeout = fitlme(lmedata, 'meas ~ 1 + (1|id)',...
            'FitMethod','REML');

        dep = depall(lmeout,trlmean);

        relsummary.group(gloc).event(eloc).dependability = dep;
        relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
        relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
        relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
        relsummary.group(gloc).event(eloc).goodn = height(goodids);
   
    case 2 %possible multiple groups but no event types to consider
        
       	eloc = 1;

        for gloc=1:ngroups

            if eloc == 1
                relsummary.group(gloc).name = gnames{gloc};
            end

            relsummary.group(gloc).event(eloc).name = enames{eloc};

            for trial=1:ntrials
                mrel(trial) = ...
                    mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,trial)); 
                llrel(trial) = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,trial),.025);
                ulrel(trial) = quantile(reliab(...
                    data.g(gloc).e(eloc).sig_u.raw,...
                    data.g(gloc).e(eloc).sig_e.raw,trial),.975);
            end

            %find the number of trials to reach cutoff based on the
            %lower limit of the confidence interval
            trlcutoff = find(llrel >= relcutoff, 1);

            relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
            relsummary.group(gloc).event(eloc).mrel = mrel(trlcutoff);
            relsummary.group(gloc).event(eloc).llrel = llrel(trlcutoff);
            relsummary.group(gloc).event(eloc).ulrel = ulrel(trlcutoff);


            %Only calculate overall dependability on the ids with
            %enough trials

            datatrls = REL.data;
            ind = strcmp(datatrls.group,gnames{gloc});
            datatrls = datatrls(ind,:);

            trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});

            ind2include = trltable.GroupCount >= trlcutoff;
            ind2exclude = trltable.GroupCount < trlcutoff;

            relsummary.group(gloc).event(eloc).eventgoodids =...
                trltable.id(ind2include);
            relsummary.group(gloc).event(eloc).eventbadids =...
                trltable.id(ind2exclude);
        end

        %find the ids that have enough trials for each event type
        
        for gloc=1:ngroups
        
            relsummary.group(gloc).goodids = ...
                relsummary.group(gloc).event(eloc).eventgoodids;
            relsummary.group(gloc).badids = ...
                relsummary.group(gloc).event(eloc).eventbadids;
            
        end

        for gloc=1:ngroups

            datatable = REL.data;
            ind = strcmp(datatable.group,gnames{gloc});
            datasubset = datatable(ind,:);

            goodids = table(relsummary.group(gloc).goodids);

            lmedata = innerjoin(datasubset, goodids,...
                'LeftKeys', 'id', 'RightKeys', 'Var1',...
                'LeftVariables', {'id' 'meas'});

            trltable = varfun(@length,lmedata,...
                'GroupingVariables',{'id'});

            trlmean = mean(trltable.GroupCount);

            %perform calculations for overall reliability
            lmeout = fitlme(lmedata, 'meas ~ 1 + (1|id)',...
                'FitMethod','REML');

            dep = depall(lmeout,trlmean);

            relsummary.group(gloc).event(eloc).dependability = dep;
            relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
            relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
            relsummary.group(gloc).event(eloc).goodn = height(goodids);
        end
        
        
    case 3 %possible event types but no groups to consider
        
        gloc = 1;
        for eloc=1:nevents
                
                if eloc == 1
                    relsummary.group(gloc).name = gnames{gloc};
                end
               
                relsummary.group(gloc).event(eloc).name = enames{eloc};
                
                for trial=1:ntrials
                    mrel(trial) = ...
                        mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial)); 
                    llrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.025);
                    ulrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.975);
                end
                
                %find the number of trials to reach cutoff based on the
                %lower limit of the confidence interval
                trlcutoff = find(llrel >= relcutoff, 1);
                
                relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                relsummary.group(gloc).event(eloc).mrel = mrel(trlcutoff);
                relsummary.group(gloc).event(eloc).llrel = llrel(trlcutoff);
                relsummary.group(gloc).event(eloc).ulrel = ulrel(trlcutoff);
                
                
                %Only calculate overall dependability on the ids with
                %enough trials
                
                datatrls = REL.data{eloc};
                
                trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                
                ind2include = trltable.GroupCount >= trlcutoff;
                ind2exclude = trltable.GroupCount < trlcutoff;
                
                relsummary.group(gloc).event(eloc).eventgoodids = trltable.id(ind2include);
                relsummary.group(gloc).event(eloc).eventbadids = trltable.id(ind2exclude);

        end
        
        tempids = {};
        badids = [];
        for eloc=1:nevents
            tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;

            if eloc > 1
                [~,ind]=setdiff(tempids{1},tempids{eloc});
                new = tempids{1};
                badids = [badids...
                    relsummary.group(gloc).event(eloc).eventbadids...
                    new(ind)];
                new(ind) = [];
                tempids{1} = new;
            end

        end

        relsummary.group(gloc).goodids = tempids{1};
        relsummary.group(gloc).badids = unique(badids);
        
        for eloc=1:nevents

                datatable = REL.data{eloc};
               
                goodids = table(relsummary.group(gloc).goodids);
                
                lmedata = innerjoin(datatable, goodids,...
                    'LeftKeys', 'id', 'RightKeys', 'Var1',...
                    'LeftVariables', {'id' 'meas'});
                
                trltable = varfun(@length,lmedata,...
                    'GroupingVariables',{'id'});
                
                trlmean = mean(trltable.GroupCount);
                
                %perform calculations for overall reliability
                lmeout = fitlme(lmedata, 'meas ~ 1 + (1|id)',...
                    'FitMethod','REML');
                
                dep = depall(lmeout,trlmean);
                
                relsummary.group(gloc).event(eloc).dependability = dep;
                relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                relsummary.group(gloc).event(eloc).goodn = height(goodids);

        end
                
    case 4 %groups and event types to consider
        
        for eloc=1:nevents
            
            for gloc=1:ngroups
                
                if eloc == 1
                    relsummary.group(gloc).name = gnames{gloc};
                end
               
                relsummary.group(gloc).event(eloc).name = enames{eloc};
                
                for trial=1:ntrials
                    mrel(trial) = ...
                        mean(reliab(data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial)); 
                    llrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.025);
                    ulrel(trial) = quantile(reliab(...
                        data.g(gloc).e(eloc).sig_u.raw,...
                        data.g(gloc).e(eloc).sig_e.raw,trial),.975);
                end
                
                %find the number of trials to reach cutoff based on the
                %lower limit of the confidence interval
                trlcutoff = find(llrel >= relcutoff, 1);
                
                relsummary.group(gloc).event(eloc).trlcutoff = trlcutoff;
                relsummary.group(gloc).event(eloc).mrel = mrel(trlcutoff);
                relsummary.group(gloc).event(eloc).llrel = llrel(trlcutoff);
                relsummary.group(gloc).event(eloc).ulrel = ulrel(trlcutoff);
                
                
                %Only calculate overall dependability on the ids with
                %enough trials
                
                datatrls = REL.data{eloc};
                ind = strcmp(datatrls.group,gnames{gloc});
                datatrls = datatrls(ind,:);
                
                trltable = varfun(@length,datatrls,'GroupingVariables',{'id'});
                
                ind2include = trltable.GroupCount >= trlcutoff;
                ind2exclude = trltable.GroupCount < trlcutoff;
                
                relsummary.group(gloc).event(eloc).eventgoodids = trltable.id(ind2include);
                relsummary.group(gloc).event(eloc).eventbadids = trltable.id(ind2exclude);

                 
            end
        end
        
        %find the ids that have enough trials for each event type
        
        for gloc=1:ngroups
            tempids = {};
            badids = [];
            for eloc=1:nevents
                tempids{end+1} = relsummary.group(gloc).event(eloc).eventgoodids;

                if eloc > 1
                    [~,ind]=setdiff(tempids{1},tempids{eloc});
                    new = tempids{1};
                    badids = [badids...
                        relsummary.group(gloc).event(eloc).eventbadids...
                        new(ind)];
                    new(ind) = [];
                    tempids{1} = new;
                end

                

            end

        relsummary.group(gloc).goodids = tempids{1};
        relsummary.group(gloc).badids = unique(badids);

        end
         
        
        for eloc=1:nevents

            for gloc=1:ngroups
                
                datatable = REL.data{eloc};
                ind = strcmp(datatable.group,gnames{gloc});
                datasubset = datatable(ind,:);
                
                goodids = table(relsummary.group(gloc).goodids);
                
                lmedata = innerjoin(datasubset, goodids,...
                    'LeftKeys', 'id', 'RightKeys', 'Var1',...
                    'LeftVariables', {'id' 'meas'});
                
                trltable = varfun(@length,lmedata,...
                    'GroupingVariables',{'id'});
                
                trlmean = mean(trltable.GroupCount);
                
                %perform calculations for overall reliability
                lmeout = fitlme(lmedata, 'meas ~ 1 + (1|id)',...
                    'FitMethod','REML');
                
                dep = depall(lmeout,trlmean);
                
                relsummary.group(gloc).event(eloc).dependability = dep;
                relsummary.group(gloc).event(eloc).trlinfo.min = min(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.max = max(trltable.GroupCount);
                relsummary.group(gloc).event(eloc).trlinfo.mean = trlmean;
                relsummary.group(gloc).event(eloc).goodn = height(goodids);
            end
        end
                
        
end %switch analysis

%reminder....
%1 - no groups or event types to consider
%2 - possible multiple groups but no event types to consider
%3 - possible event types but no groups to consider
%4 - possible groups and event types to consider

%Create table to display both sets of data
label = {};
trlcutoff = {};
dep = {};
overalldep = {};
mintrl = {};
maxtrl = {};
meantrl = {};
goodn = {};
badn = {};

%put data together in tables to display

for gloc=1:ngroups
    for eloc=1:nevents
        switch analysis
            case 1
                label{end+1} = 'Measurement';
            case 2
                label{end+1} = gnames{gloc};
            case 3
                label{end+1} = enames{eloc};
            case 4
                label{end+1} = [gnames{gloc} ' - ' enames{eloc}];
        end
        trlcutoff{end+1} = relsummary.group(gloc).event(eloc).trlcutoff;
        dep{end+1} = sprintf('        %0.2f CI[ %0.2f, %0.2f]',...
            relsummary.group(gloc).event(eloc).mrel,...
            relsummary.group(gloc).event(eloc).llrel,...
            relsummary.group(gloc).event(eloc).ulrel);
        
        overalldep{end+1} = ...
            str2num(sprintf('%0.2f',...
            round(relsummary.group(gloc).event(eloc).dependability,2)));
        mintrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.min;
        maxtrl{end+1} = relsummary.group(gloc).event(eloc).trlinfo.max;
        meantrl{end+1} = ...
            round(relsummary.group(gloc).event(eloc).trlinfo.mean,2);
        goodn{end+1} = relsummary.group(gloc).event(eloc).goodn;
        badn{end+1} = length(relsummary.group(gloc).badids);
        
    end 
end

inctrltable = table(label',trlcutoff',dep');

overalltable = table(label',goodn',badn',overalldep',meantrl',...
    mintrl',maxtrl');


%define parameters for figure position
figwidth = 429;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

ra_inctrl= figure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Results of Increasing Trials on Dependability',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(ra_inctrl,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Dependability Analyses',...
    'Position',[0 row figwidth 25]);          

%Start a table
t = uitable('Parent',ra_inctrl,'Position',...
    [25 (rowspace*2) figwidth-50 figheight-(rowspace*5)],...
    'Data',table2cell(inctrltable));
set(t,'ColumnName',{'Label' 'Trial Cutoff' 'Dependability'});
set(t,'ColumnWidth',{150 'auto' 150});
set(t,'RowName',[]);
set(t,'FontSize',12);

%define parameters for figure position
figwidth = 650;
figheight = 400;

%define space between rows and first row location
rowspace = 25;
row = figheight - rowspace*2;

ra_overall= figure('unit','pix',...
  'position',[400 400 figwidth figheight],...
  'menub','no',...
  'name','Dependability Analyses Including All Trials',...
  'numbertitle','off',...
  'resize','off');

%Print the name of the loaded dataset
uicontrol(ra_overall,'Style','text','fontsize',16,...
    'HorizontalAlignment','center',...
    'String','Overall Dependability',...
    'Position',[0 row figwidth 25]);          

%Start a table
t = uitable('Parent',ra_overall,'Position',...
    [25 (rowspace*2) figwidth-50 figheight-(rowspace*5)],...
    'Data',table2cell(overalltable));
set(t,'ColumnName',{'Label' 'n Included' 'n Excluded' ...
    'Dependability' 'Mean # of Trials' 'Min # of Trials'...
    'Max # of Trials'});
set(t,'RowName',[]);
set(t,'FontSize',12);

RELout = REL;
RELout.relsummary = relsummary;

end

function depout = reliab(var_u, var_e, obs)

depout = var_u.^2 ./ (var_u.^2 + (var_e.^2./obs));

end

function iccout = icc(var_u,var_e)

iccout = var_u.^2 ./ (var_u.^2 + var_e.^2);

end

function dep = depall(lme,num)

[obsvar,resvar] = covarianceParameters(lme);

dep = cell2mat(obsvar)/(cell2mat(obsvar)+(resvar/num));

end