function ra_relfigures(varargin)
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

%TODOS
%ensure that everything works with 1 group or 1 event or no group or event
%add the overall level of dependability after implementing given cutoffs
%store dependability estimates in REL structure array

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
% data.g.glabel = {};
% data.g.elabel = {};
% data.g.mu = [];
% data.g.mu.raw = [];
% data.g.mu.ll = [];
% data.g.mu.ul = [];
% data.g.sig_u = [];
% data.g.sig_u.raw = [];
% data.g.sig_u.ll = [];
% data.g.sig_u.ul = [];
% data.g.sig_e = [];
% data.g.sig_e.raw = [];
% data.g.sig_e.ll =[];
% data.g.sig_e.ul = [];

if strcmpi('REL.groups','none')
    ngroups = 1;
    gnames = REL.groups;
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

if strcmpi('REL.events','none')
    nevents = 1;
    enames = REL.events;
else
    nevents = length(REL.events);
    enames = REL.events(:);
end

for i=1:length(REL.out.labels)
    
    lblstr = strsplit(REL.out.labels{i},'_');
    
    eloc = find(ismember(enames,lblstr(1)));
    gloc = find(ismember(gnames,lblstr(2)));
    
    if isempty(eloc); eloc = 1; end;
    if isempty(gloc); gloc = 1; end;
        
    data(gloc).g(eloc).label = REL.out.labels(i);
    data(gloc).g(eloc).mu.raw = REL.out.mu(:,i);
    data(gloc).g(eloc).sig_u.raw = REL.out.sig_u(:,i);
    data(gloc).g(eloc).sig_e.raw = REL.out.sig_e(:,i);
    data(gloc).g(eloc).elabel = enames(eloc);
    data(gloc).g(eloc).glabel = gnames(gloc);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fitdata = table(mu, sig_u, sig_e, lp);
% 
% mean_icc = mean(icc(sig_u(:), sig_e(:)));
% ll_icc = quantile(icc(sig_u(:), sig_e(:)),.025);
% ul_icc = quantile(icc(sig_u(:), sig_e(:)),.975);
% 
% mean_sig_u = mean(sig_u);
% mean_sig_e = mean(sig_e);
% mean_mu = mean(mu);
% 
% ll_sig_u = quantile(sig_u,.025);
% ul_sig_u = quantile(sig_u,.975);
% ll_sig_e = quantile(sig_e,.025);
% ul_sig_e = quantile(sig_e,.975);
% ll_mu = quantile(mu,.025);
% ul_mu = quantile(mu,.975);


%flexible number of trials to display (warn about extrapolating too far!)
ntrials = 50;
x = 1:ntrials;

mrel = zeros(0,ntrials);
llrel = zeros(0,ntrials);
ulrel = zeros(0,ntrials);

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
                mean(reliab(data(gloc).g(eloc).sig_u.raw,...
                data(gloc).g(eloc).sig_e.raw,trial)); 
        end
    end
    subplot(yplots,xplots,eloc);
    plot(x,mrel);
    axis([0 ntrials 0 1]);
    set(gca,'fontsize',16);
    title(enames{eloc},'FontSize',20);
    ylabel('Dependability','FontSize',fsize);
    xlabel('Number of Observations','FontSize',fsize);
    hline = refline(0,.7);
    set(hline,'Color','b','LineStyle',':');
    leg = legend(gnames{:},'Location','southeast');
    set(leg,'FontSize',fsize);
end

%display in the command window the number of trials that should be included
%to achieve a given level (.70) of dependability

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

relsummary.relcutoff = relcutoff;

switch analysis
    case 1 %no groups or event types to consider
    case 2 %possible multiple groups but no event types to consider
    case 3 %possible event types but no groups to consider
        
        
        
    case 4 %groups and event types to consider
        
        for eloc=1:nevents
            
            for gloc=1:ngroups
                
                if eloc == 1
                    relsummary.group(gloc).name = gnames{gloc};
                end
               
                relsummary.group(gloc).event(eloc).name = enames{eloc};
                
                for trial=1:ntrials
                    mrel(trial) = ...
                        mean(reliab(data(gloc).g(eloc).sig_u.raw,...
                        data(gloc).g(eloc).sig_e.raw,trial)); 
                    llrel(trial) = quantile(reliab(...
                        data(gloc).g(eloc).sig_u.raw,...
                        data(gloc).g(eloc).sig_e.raw,trial),.025);
                    ulrel(trial) = quantile(reliab(...
                        data(gloc).g(eloc).sig_u.raw,...
                        data(gloc).g(eloc).sig_e.raw,trial),.975);
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
                
                relsummary.group(gloc).event(eloc).eventgoodids = trltable.id(ind2include);
                 
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
                    badids = [badids new(ind)];
                    new(ind) = [];
                    tempids{1} = new;
                end

                

            end

        relsummary.group(gloc).goodids = tempids{1};
        relsummary.group(gloc).badids = badids;

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
fprintf('\n\n\n');

if ~isempty(trlcutoff)
                    fprintf(...
                        '\n%s %s: %d trials, dependability %0.2f CI[%0.2f, %0.2f]',...
                        enames{eloc},gnames{gloc},...
                        trlcutoff,mrel(trlcutoff),...
                        llrel(trlcutoff),ulrel(trlcutoff));
                elseif isempty(trlcutoff)
                    fprintf(...
                        '\n%s %s: Level of dependability not obtained',...
                        enames{eloc},gnames{gloc});
                end



%%%%%%load file from R to make sure that I am using the correct variance
%%%%%%estimates

%perform calculations separately for each group and event
lmeout = fitlme(REL.data{1,1}, 'meas ~ 1 + (1|id)','FitMethod','REML');

%can try converting the data to a table to filter data by number of trials
%grpstats(CARel,'Participant') will give number of rows per participant

[obsvar,resvar] = covarianceParameters(lmeout);



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