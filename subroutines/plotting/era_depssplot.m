function depplot = era_depssplot(varargin)
%Plot the dependability as a function of the number of trials included in
%an average (stratified by group and condition, if applicable)
%
%depplot = era_depvtrialsplot('era_data',era_data,'trials',[1 50],...
%   'depline',plotdepline,'depcutoff',depcutoff);
%
%Last Modified 1/19/17
%
%Inputs
% era_data - ERA Toolbox data structure array. Variance components should
%  be included.
% trials - range of trials to plot
% depcutoff - dependality threshold for deeming data reliable (used for
%  reference line in plot)
%
%Optional Input
% depline - dependabilty estimate to plot (default: 2)
%  1 - lower limit
%  2 - point estimate
%  3 - upper limit
% CI - confidence interval width. Decimal from 0 to 1. (default: .95)
%
%
%Outputs
% depplot - figure handle for the plot
% figure that displays the relationship between dependability and the
%  number of trials retained for averaging stratified by group and
%  condition

% Copyright (C) 2016-2020 Peter E. Clayson
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
% by Peter Clayson (7/24/16)
% peter.clayson@gmail.com
%
%
%9/16/16 PC
% remove box around legend
%
%1/19/17 PC
% updated copyright

%somersault through inputs
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help era_depssplot for more information about inputs'));
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
            'See help era_depssplot for more information \n'));
    end
    
end


%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
ciperc = era_data.relsummary.ciperc;
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

%see how many subplots are needed
if (nevents * ngroups) > 2
    xplots = ceil(sqrt((nevents * ngroups)));
    yplots = ceil(sqrt((nevents * ngroups)));
elseif (nevents * ngroups) == 2
    xplots = 1;
    yplots = 2;
elseif (nevents * ngroups) == 1
    xplots = 1;
    yplots = 1;
end

%create the figure
depplot = figure;
set(gcf,'NumberTitle','Off');
depplot.Position = [125 630 900 450];
fsize = 16;

plottitle = 'Subject-Level Dependability Estimates with Credible Intervals';

%extract the data and create the subplots for depplot
for eloc=1:nevents
    for gloc=1:ngroups
        gtable = era_data.relsummary.group(gloc).event(eloc).ssrel_table;
        gtable = sortrows(gtable,'dep_pt');
        
        plotssrel.subj = 1:height(gtable);
        plotssrel.dep_pt = gtable.dep_pt;
        plotssrel.dep_ll = gtable.dep_pt - gtable.dep_ll;
        plotssrel.dep_ul = gtable.dep_ul - gtable.dep_pt;
        
        gro_est = gtable.bp_var(1) / ...
            (gtable.bp_var(1) + (gtable.pop_errvar(1)/mean(gtable.trls)));
        
        miss = 0;
        miss_ids = [];
        
        if any(gtable.dep_ul < gro_est) ||... 
                any(gtable.dep_ll > gro_est)
            
            miss = 1;
            miss_ids = find(gtable.dep_ul < gro_est);
            miss_ids = [miss_ids find(gtable.dep_ll > gro_est)];
            
            good_ids = 1:height(gtable);
            good_ids(miss_ids) = [];
        end
        
        
        set(gcf,'Name',plottitle);
        subplot(yplots,xplots,eloc);
        if miss
            hold on
            errorbar(plotssrel.subj(miss_ids),...
                plotssrel.dep_pt(miss_ids),...
                plotssrel.dep_ll(miss_ids),...
                plotssrel.dep_ul(miss_ids),...
                'Marker','.',...
                'MarkerSize',25,...
                'LineWidth',1.5,...
                'LineStyle','none',...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'Color','r');
            
            errorbar(plotssrel.subj(good_ids),...
                plotssrel.dep_pt(good_ids),...
                plotssrel.dep_ll(good_ids),...
                plotssrel.dep_ul(good_ids),...
                'Marker','.',...
                'MarkerSize',25,...
                'LineWidth',1.5,...
                'LineStyle','none',...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','b',...
                'Color','b');
            
            hold off
        else
            errorbar(plotssrel.subj,...
            plotssrel.dep_pt,...
            plotssrel.dep_ll,...
            plotssrel.dep_ul,...
            'Marker','.',...
            'MarkerSize',25,...
            'LineWidth',1.5,...
            'LineStyle','none');
            
        end
        axis([0 max(plotssrel(gloc).subj)+1 0 1]);
        set(gca,'fontsize',16);
        
        if (~strcmpi(enames{eloc},'none') &&... 
                ~strcmpi(gnames{gloc},'none')) &&... 
                nevents > 1 && ngroups > 1
            title([gnames{gloc} ': ' enames{eloc}],'FontSize',20);
        elseif (strcmpi(enames{eloc},'none') &&... 
                ~strcmpi(gnames{gloc},'none')) && nevents > 1
            title(enames{eloc},'FontSize',20);
        elseif (~strcmpi(enames{eloc},'none') &&... 
                strcmpi(gnames{gloc},'none')) && ngroups > 1
            title(gnames{gloc},'FontSize',20);
        else 
            title('');
        end
        
        ylabel('Dependability','FontSize',fsize);
        xlabel('Participants Ordered by Dependability Estimate',...
            'FontSize',fsize);
        hline = refline(0,gro_est);
        set(hline,...
            'Color','k',...
            'LineStyle','-',...
            'LineWidth',1.5);
    end
end


end