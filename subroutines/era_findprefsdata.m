function [era_prefs,era_data] = era_findprefsdata(varin)
%
%Somersault through inputs to find era_prefs and era_data
%
%[era_prefs,era_data] = era_findprefsdata(varargin)
%
%Last Updated 7/21/16 
%
%Input
% varin - varargin from guis
%
%output
% era_prefs - ERA Toolbox structure array containing preferences
% era_data - ERA Toolbox structure array containing data
%

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
% by Peter Clayson (7/21/16)
% peter.clayson@gmail.com
%

%Somersault through varargin inputs to check for which inputs were
%defined and store those values. 
%check if era_prefs has been defined
ind = find(strcmp('era_prefs',varin),1);
if ~isempty(ind)
    era_prefs = varin{ind+1}; 
else
    era_prefs = '';
end

%check if era_data has been defined
ind = find(strcmp('era_data',varin),1);
if ~isempty(ind)
    era_data = varin{ind+1}; 
else
    era_data = '';
end

end

