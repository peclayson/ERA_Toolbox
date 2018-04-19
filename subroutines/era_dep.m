function [ll,pt,ul] = era_dep(varargin)
%Calculate dependability using variance components from CmdStan
%
%[ll,pt,ul] = era_dep('bp',var_u,'wp',var_e,'obs',[1 50],'CI',.95)
%
%Last Modified 1/19/17
%
%Inputs
% bp - between-person variance components from CmdStan (var_u)
% wp - with-person variance components from CmdStan (var_e)
% obs - number of observations. Can either be 1 number or range [min max]
%  E.g., The overall dependability estimate will use one number (central
%   tendency). The plot for the number of trials v dependability will use a
%   range of numbers.
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Outputs
% ll - lower limit of the credible interval specified by CI for
%  dependability
% pt - the point estimate for dependability
% ul - upper limit of the credible interval specified by CI for
%  dependability

% Copyright (C) 2016-2018 Peter E. Clayson
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
% by Peter Clayson (7/23/16)
% peter.clayson@gmail.com
%
%1/19/17 PC
% updated copyright
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
        'See help era_dep for more information about inputs'));
    end
    
    %check if bp was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1}; 
    else 
        error('varargin:bp',... %Error code and associated error
        strcat('WARNING: Between-person variance not specified \n\n',... 
        'Please input bp. See help era_dep for more information \n'));
    end
    
    %check if wp was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('wp',varargin),1);
    if ~isempty(ind)
        wp = varargin{ind+1}; 
    else 
        error('varargin:wp',... %Error code and associated error
        strcat('WARNING: Within-person variance not specified \n\n',... 
        'Please input wp. See help era_dep for more information \n'));
    end
    
    %check if obs was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('obs',varargin),1);
    if ~isempty(ind)
        obs = varargin{ind+1}; 
    else 
        error('varargin:obs',... %Error code and associated error
        strcat('WARNING: Number of observations not specified \n\n',... 
        'Please input obs. See help era_dep for more information \n'));
    end
    
    %check if CI was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help era_dep for more information \n'));
        end
    else 
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',... 
            'Please input CI. See help era_dep for more information \n'));
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

%calculate edges for CI
ciedge = (1-ciperc)/2;

%calculate the dependability for a single number of trials if obs only
%contains one number
if length(obs) == 1
    
    %dependability formula
    dep = bp.^2 ./ (bp.^2 + (wp.^2./obs));
    
    %pull the lower limit, point estimate, and upper limit
    ll = quantile(dep,ciedge);
    pt = mean(dep);
    ul = quantile(dep,1-ciedge);
    
elseif length(obs) == 2

    %if obs contains two numbers calculate the dependabiltiy for the range
    %specified
    
    %create empty matrices for storing values
    ll = zeros(0,obs(2)-obs(1)+1);
    pt = zeros(0,obs(2)-obs(1)+1);
    ul = zeros(0,obs(2)-obs(1)+1);
    
    for i = obs(1):obs(2)
        %dependability formula
        dep = bp.^2 ./ (bp.^2 + (wp.^2./i));

        %pull the lower limit, point estimate, and upper limit
        ll(i) = quantile(dep,ciedge);
        pt(i) = mean(dep);
        ul(i) = quantile(dep,1-ciedge);  
        
    end %for i = obs(1):obs(2)
    
    
end %if length(obs) == 1

end