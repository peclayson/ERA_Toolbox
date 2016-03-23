function result = era_checkconvergence(fit)
%
%Check r_hat values for convergence from a StanFit model
%
%era_checkconvergence(fit)
%
%Lasted Updated 3/23/16
%
%Required Input:
% fit - Stanfit model output
%
%Outputs:
% result - 1, chains converge; 0, chains failed to converge

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
% by Peter Clayson (3/23/16)
% peter.clayson@gmail.com

%pull summary using the print function
%there's no other way to get at the r_hat values, so output will also be
%printed in the command window
output = print(fit);

%define cell array that holds the r_hats
allrhats = cell(0,2);

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
        rhat = str2num(pullrow{end});
        allrhats{end+1,1} = label;
        allrhats{end,2} = rhat;
    end
end

%see if any of the rhats didn't equal 1 (i.e., did not converge)
indbad = find([allrhats{:,2}] ~= 1);

%create a result output that we can use to store convergence values
result = struct;

%specify whethere convergenece between chains was reached
if ~isempty(indbad)
    result.converged = 0;
else
    result.converged = 1;
end

%store the convergence values
result.rhats = allrhats;

%if the results didn't converge recommend that the user re-run
if result.converged == 0
    
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
    
    str = {'Chains did not converge';...
        'Would you like to rerun with more iterations?'};
    
    %Write text
    uicontrol(era_gui,'Style','text','fontsize',16,...
        'HorizontalAlignment','center',...
        'String',str,...
        'Position',[0 row figwidth 50]);          

    %Create a button that will take the user to the gui for setting the inputs
    %to process data
    uicontrol(era_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','<html><center>Do Not<br>Rerun',...
        'Position', [figwidth/8 25 figwidth/3 75],...
        'Callback',{@era_conv_startview,era_gui}); 

    %Create button that will take the user to the gui for setting the inputs
    %for viewing the data
    uicontrol(era_gui,'Style','push','fontsize',14,...
        'HorizontalAlignment','center',...
        'String','Rerun',...
        'Position', [5*figwidth/8 25 figwidth/3 75],...
        'Callback',{@era_conv_recompute}); 

    %display gui
    set(era_gui,'Visible','on');
end
    
   
end

function era_conv_startview(varargin)

close(varargin{end});

assign('era_computerel','converged',3);

end


function era_conv_recompute(varargin)

close(varargin{end});

assign('era_computerel','iter',3);

end
