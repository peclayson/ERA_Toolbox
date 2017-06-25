function depvers = era_dependentsversions
%
%Specifies the versions used for the dependents of the ERA Toolbox
%
%Last Updated 6/24/17
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

% Copyright (C) 2016-2017 Peter E. Clayson
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

depvers = struct;
depvers.cmdstan = '2.16.0';
depvers.matlabstan = '2.15.1.0';
depvers.matlabprocessmanager = '0.5.1';

end

