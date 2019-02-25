function era_checktraceplots(REL,varargin)
%View trace plots of Stan parameters and ask the user whether Stan should
% rerun or save the outputs
%
%Last Updated 8/23/17
%
%To identify whether the estimates are stable, the should look like "fat,
%hairy caterpillars". See User Manual for more information.
%
%Input:
% REL - output from stan after running data through era_computerel
%
%Optional Input:
% askuser - displays a gui asking the user whether to rerun the model. 0 -
%  no (default); 1 - yes
%
%Output:
% No output generated. The user's choice will be checked in era_startproc
%

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
% by Peter Clayson (8/17/17)
% peter.clayson@gmail.com
%
%8/22/17 PC
% changes after running the checktraceplots function from the gui
% new input to specify whether the user should be asked about rerunning the
%  model
%
%8/23/17 PC
% after continued testing I realized that it was collapsing across group
%  and event for each parameter - fixed
% additional to fix to counting groups and events

%check for era_prefs
if nargin == 0 || isempty(REL)
    error('REL:notfound',...
        strcat('ERROR: Input is incomplete \n\n',...
        'REL from the era_computerel function is needed\n',...
        'See help for era_checktraceplots more information'));
end

if ~isempty(varargin)
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Optional inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_relfigures for more information on optional inputs'));
    end
    
    %check whether askuser was provided
    ind = find(strcmp('askuser',varargin),1);
    if ~isempty(ind)
        askuser = varargin{ind+1}; 
    else 
        askuser = 0;
    end
end

%check whether any groups exist
if strcmpi(REL.groups,'none')
    ngroups = 1;
    %gnames = cellstr(era_data.rel.groups);
    gnames ={''};
else
    ngroups = length(REL.groups);
    gnames = REL.groups(:);
end

%check whether any events exist
if strcmpi(REL.events,'none')
    nevents = 1;
    %enames = cellstr(REL.events);
    enames = {''};
else
    nevents = length(REL.events);
    enames = REL.events(:);
end

%create an x-axis for the number of observations
x = 1:REL.niter/2;

%create an empty array for storing information into
tpdata = zeros(REL.niter/2,3);

%define how many subplots are needed
xplots = 3;
yplots = nevents * ngroups;

%create the figure
tplots = figure;
set(gcf,'NumberTitle','Off');
tplots.Position = [125 1100 1000 1100];
fsize = 16;
set(gcf,'Name', 'Trace Plots of Stan Parameters');

trackwhichplot = 1;
countgroupevent = 1;
nrow = 1;

%Create the trace plots for each parameter
for nplot=1:(yplots * 3)
    subplot(yplots,xplots,nplot); 
    axis([0 REL.niter/2 0 1]);
    
    %keep track of which plot is being plotted so it ends up in the
    %appropriate column with the appropriate labels
    if trackwhichplot == 1
        
        for jj = 1:REL.nchains
            tpdata(:,jj) = REL.out.mu(((jj-1)*(REL.niter/2))+1:...
                jj*(REL.niter/2),nrow);
        end
        
        plot(x,tpdata);
        
        if contains(REL.out.labels,'_;_')
            str = strsplit(REL.out.labels{countgroupevent},'_;_');
            str = [str{1} ' - ' str{2}];
        else
            str = REL.out.labels;
        end

        ylabel(str,'FontSize',fsize);

        printtitles(yplots,nplot,fsize);
        
        trackwhichplot = trackwhichplot + 1;

    elseif trackwhichplot == 2

        for jj = 1:REL.nchains
            tpdata(:,jj) = REL.out.sig_u(((jj-1)*(REL.niter/2))+1:...
                jj*(REL.niter/2),nrow);
        end

        plot(x,tpdata);
        
        printtitles(yplots,nplot,fsize);
        
        trackwhichplot = trackwhichplot + 1;

    elseif trackwhichplot == 3

        for jj = 1:REL.nchains
            tpdata(:,jj) = REL.out.sig_e(((jj-1)*(REL.niter/2))+1:...
                jj*(REL.niter/2),nrow);
        end

        plot(x,tpdata);
        
        printtitles(yplots,nplot,fsize);
        
        trackwhichplot = 1;
        nrow = nrow + 1;
        countgroupevent = countgroupevent + 1;
    end

end

%tag the gui so it can be closed in  era_computerelwrap
tplots.Tag = 'tplots';

if askuser == 1
    %define parameters for figure position
    figwidth = 500;
    figheight = 200;

    %define space between rows and first row location
    rowspace = 25;
    row = figheight - rowspace*2.5;

    %initialize gui
    era_gui= figure('unit','pix','Visible','off',...
      'position',[400 400 figwidth figheight],...
      'menub','no',...
      'numbertitle','off',...
      'resize','off');

    movegui(era_gui,'center');

    str = {'Would you like to rerun with more iterations?'};

    %Write text
    uicontrol(era_gui,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',str,...
        'Position',[0 row figwidth 50]);          

    %Create a button that will not rerun the model
    uicontrol(era_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','<html><center>Do Not<br>Rerun',...
        'Position', [figwidth/8 25 figwidth/3 75],...
        'Callback',{@era_setrerun,0,era_gui}); 

    %Create button that will rerun the model with more iterations
    uicontrol(era_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Rerun',...
        'Position', [5*figwidth/8 25 figwidth/3 75],...
        'Callback',{@era_setrerun,1,era_gui}); 

    %display gui
    set(era_gui,'Visible','on');

    %tag the gui
    era_gui.Tag = 'era_gui';

    %wait so the user can respond
    uiwait(era_gui);
end

end

function era_setrerun(varargin)
%to rerun or not to rerun the model
uiresume;
era_gui = findobj('Tag','era_gui');
guidata(era_gui,varargin{3});

end

function printtitles(yplots,nplot,fsize)
%functiont to print axis titles
if nplot == 1
    title('Mean','FontSize',fsize);
elseif nplot == 2
    title('Between-Person Variance','FontSize',fsize);
elseif nplot == 3
    title('Within-Person Variance','FontSize',fsize);
end
    
if (yplots * 3 - 2) <= nplot
    xlabel('Number of Samples','FontSize',fsize);
end
    
end




