function ra_relfigures(varargin)
%Creates figures depicting dependability estimates and displays information
% about optimal cutoffs and overall dependability
%
%ra_relfigures('data',REL)
%
%Required Inputs:
% data - structure array containing results of dependability analyses using
%  cmdstan. See ra_computerel for more information.
%
%Output:
% One figure is plotted that displays the dependability of measurements as
%  the number of trials increases. If there is separate information for
%  groups, events, or groups and events, the dependability estimates will
%  be displayed accordingtly.
%
%History 
% by Peter Clayson (12/23/15)
% peter.clayson@gmail.com

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
    hline.Color = 'b';
    hline.LineStyle = ':';
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

fprintf('\n\n\nIn order to achieve a .70 level of dependability');
fprintf('\nthe following cutoffs should be used\n');

switch analysis
    case 1 %no groups or event types to consider
    case 2 %possible multiple groups but no event types to consider
    case 3 %possible event types but no groups to consider
    case 4 %groups and event types to consider
        for eloc=1:nevents
            for gloc=1:ngroups
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
                trlcutoff = find(mrel >= .70, 1);
                fprintf('\n%s %s: %d trials, dependability CI[%0.2f, %0.2f]',...
                    enames{eloc},gnames{gloc},...
                    trlcutoff,llrel(trlcutoff),ulrel(trlcutoff));
            end
        end
end
fprintf('\n\n\n');


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
