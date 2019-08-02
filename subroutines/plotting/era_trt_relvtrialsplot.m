function relplot = era_trt_relvtrialsplot(varargin)
%Plot the reliability as a function of the number of trials included in
%an average (stratified by group and condition, if applicable)
%
%relplot = era_trt_relvtrialsplot('era_data',era_data,'trials',[1 50],...
%   'relline',plotrelline,'relcutoff',relcutoff);
%
%Last Modified 8/2/19
%
%Inputs
% era_data - ERA Toolbox data structure array. Variance components should
%  be included.
% trials - range of trials to plot
% relcutoff - reliability threshold for deeming data reliable (used for
%  reference line in plot)
%
%Optional Input
% relline - reliability estimate to plot (default: 2)
%  1 - lower limit
%  2 - point estimate
%  3 - upper limit
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% relplot - figure handle for the plot
% figure that displays the relationship between dependability and the
%  number of trials retained for averaging stratified by group and
%  condition

% Copyright (C) 2016-2019 Peter E. Clayson
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
% by Peter Clayson (8/2/19)
% peter.clayson@gmail.com
%


%somersault through inputs
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_trt_relvtrialsplot for more information about inputs'));
    end
    
    %check if era_data was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('era_data',varargin),1);
    if ~isempty(ind)
        era_data = varargin{ind+1}; 
    else 
        error('varargin:era_data',... %Error code and associated error
            strcat('WARNING: era_data not specified \n\n',... 
            'Please input era_data (ERA Toolbox data structure array).\n',...
            'See help era_trt_relvtrialsplot for more information \n'));
    end
    
    %check if trials was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('trials',varargin),1);
    if ~isempty(ind)
        trials = varargin{ind+1}; 
    else 
        error('varargin:trials',... %Error code and associated error
            strcat('WARNING: Range of trials to plot not specified \n\n',... 
            'Please input trials: [min max]\n',...
            'See help era_trt_relvtrialsplot for more information\n'));
    end
    
    %check if relcutoff was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('relcutoff',varargin),1);
    if ~isempty(ind)
        relcutoff = varargin{ind+1}; 
    else 
        error('varargin:relcutoff',... %Error code and associated error
            strcat('WARNING: Reliability threshold not specified \n\n',... 
            'Please input relcutoff\n',...
            'See help era_trt_relvtrialsplot for more information\n'));
    end
    
    %check if relline was specified. 
    %If it is not found, set as default: 2.
    ind = find(strcmpi('relline',varargin),1);
    if ~isempty(ind)
        relline = varargin{ind+1}; 
        %make sure relline is 1, 2, or 3
        if ~any(relline==[1 2 3])
            error('varargin:relline',... %Error code and associated error
                strcat('WARNING: Reliability estimate to plot not specified \n\n',... 
                'Please enter a valid input for relline\n',...
                '1 - lower limit of credible interval\n',...
                '2 - point estimate of credible interval\n',...
                '3 - upper limit of credible interval\n',...
                'See help era_trt_relvtrialsplot for more information \n'));
        end
    else 
        relline = 2; %default is point estimate
    end
    
    %check if CI was specified. 
    %If it is not found, set default as 95%
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help era_trt_relvtrialsplot for more information \n'));
        end
    else 
        ciperc = .95; %default: 95%
    end
    
end

%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
if ciperc < .5 
    str = sprintf(' %2.0f%%',100*ciperc);
    warning('ci:width',... %Warning code and associated warning
        strcat('WARNING: Size of credible interval is small \n\n',...
        'User specified a credible interval width of ',str,'%\n',...
        'If this was intended, ignore this warning.\n'));
end

%check whether any groups exist
if strcmpi(era_data.rel.groups,'none')
    ngroups = 1;
    %gnames = cellstr(era_data.rel.groups);
    gnames ={''};
else
    ngroups = length(era_data.rel.groups);
    gnames = era_data.rel.groups(:);
end

%check whether any events exist
if strcmpi(era_data.rel.events,'none')
    nevents = 1;
    %enames = cellstr(era_data.rel.events);
    enames = {''};
else
    nevents = length(era_data.rel.events);
    enames = era_data.rel.events(:);
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

%create an x-axis for the number of observations
x = trials(1):trials(2);

%create an empty array for storing information into
plotrel = zeros(trials(2)-trials(1),0);

%see how many subplots are needed
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

%create the figure
relplot = figure;
set(gcf,'NumberTitle','Off');
relplot.Position = [125 630 900 450];
fsize = 16;

%extract the data and create the subplots for depplot
for eloc=1:nevents
    for gloc=1:ngroups
        switch relline
            case 1 %lower limit
                [plotrel(x,gloc),~,~] = era_rel_trt(...
                    'gcoeff',era_data.relsummary.gcoeff,...
                    'reltype',era_data.relsummary.reltype,...
                    'bp',era_data.relsummary.data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',era_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',era_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',era_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',era_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',era_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',era_data.relsummary.data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[trials(1) trials(2)],...
                    'CI',ciperc);
                plottitle = 'Lower Limit of 95% Credible Interval';
            case 2 %point estimate
                [~,plotrel(x,gloc),~] = era_rel_trt(...
                    'gcoeff',era_data.relsummary.gcoeff,...
                    'reltype',era_data.relsummary.reltype,...
                    'bp',era_data.relsummary.data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',era_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',era_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',era_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',era_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',era_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',era_data.relsummary.data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[trials(1) trials(2)],...
                    'CI',ciperc);
                plottitle = 'Point Estimate';
            case 3 %upper limit
                [~,~,plotrel(x,gloc)] = era_rel_trt(...
                    'gcoeff',era_data.relsummary.gcoeff,...
                    'reltype',era_data.relsummary.reltype,...
                    'bp',era_data.relsummary.data.g(gloc).e(eloc).sig_id.raw,...
                    'bo',era_data.relsummary.data.g(gloc).e(eloc).sig_occ.raw,...
                    'bt',era_data.relsummary.data.g(gloc).e(eloc).sig_trl.raw,...
                    'txp',era_data.relsummary.data.g(gloc).e(eloc).sig_trlxid.raw,...
                    'oxp',era_data.relsummary.data.g(gloc).e(eloc).sig_occxid.raw,...
                    'txo',era_data.relsummary.data.g(gloc).e(eloc).sig_trlxocc.raw,...
                    'err',era_data.relsummary.data.g(gloc).e(eloc).sig_err.raw,...
                    'obs',[trials(1) trials(2)],...
                    'CI',ciperc);
                plottitle = 'Upper Limit of 95% Credible Interval';
        end
    end
    
    switch era_data.relsummary.reltype_name
        case 'ic'
            pref1 = 'Coefficient of Equivalence, ';
        case 'trt'
            pref1 = 'Coefficient of Stability, ';
    end
    
    switch era_data.relsummary.gcoeff_name
        case 'dep'
            pref = 'Dependability v Number of Trials: ';
        case 'gen'
            pref = 'Generalizability v Number of Trials: ';
    end
    
    set(gcf,'Name',[pref1 pref plottitle]);
    subplot(yplots,xplots,eloc); 
    h = plot(x,plotrel);
    
    axis([0 trials(2) 0 trials(1)]);
    set(gca,'fontsize',16);
    
    if ~strcmpi(enames{eloc},'none')
        title(enames{eloc},'FontSize',20);
    end
    
    switch era_data.relsummary.gcoeff_name
        case 'dep'
            ylabel('Dependability','FontSize',fsize);
        case 'gen'
            ylabel('Generalizability','FontSize',fsize);
    end
    
    xlabel('Number of Observations','FontSize',fsize);
    hline = refline(0,relcutoff);
    set(hline,'Color','b','LineStyle',':');
    if analysis ~= 1 && analysis ~= 3
        leg = legend(gnames{:},'Location','southeast');
        set(leg,'FontSize',fsize);
        
        legend('boxoff');
    end

end

end