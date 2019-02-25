function depvers = era_dependentsversions
%Specifies the versions used for the dependents of the ERA Toolbox
%
%Last Updated 9/7/17
%

%This function just sets the version numbers for CmdStan, MatlabStan, and
% MatlabProcessManager
%
%Input
% No inputs required in the command line
%
%Output
% depvers - versions for each dependent
%   cmdstan - version of CmdStan
%   matlabstan - version of MatlabStan
%   matlabprocessmanager - version of MatlabProcessManager
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
% by Peter Clayson (6/24/17)
% peter.clayson@gmail.com
%
%9/7/17
% Update cmdstan version. I will not record here each time a version is
%  updated, but I will record if I make changes to the script itself.

depvers = struct;
depvers.cmdstan = '2.18.0';
depvers.matlabstan = '2.15.1.0';
depvers.matlabprocessmanager = '0.5.1';

end

