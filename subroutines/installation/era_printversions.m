function era_printversions
%
%Print the versions of Matlab, ERA Toolbox, and dependents to command
%window
%
%Last Updated 8/15/17
%

%This function just will examine whether the dependents are up to date. It
%currently does not work for MatlabStan because the version info is not
%contained in any of the downloaded files.
%
%Input
% No inputs required in the command line
%
%Output
% The versions are printed in the command window
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
% by Peter Clayson (8/15/17)
% peter.clayson@gmail.com
%

%create an empty array for storing the information about versions
versions = struct;
versions.matlab = [];
versions.ERA = [];
versions.cmdstan = [];
versions.matlabstan = [];
versions.matlabprocessmanager = [];

%store verion of Matlab
versions.matlab = version;

%store version of the ERA Toolbox
versions.ERA = era_defineversion;

%check the version of cmdstan
%find the file that has the information about the version
cs_filepath = which('cmdstan-guide.tex');

%load the file and pull the version information
cs_file = fileread(cs_filepath);
C = strsplit(cs_file,'\');
str = C{5};
c_beg = strfind(str,'{');
versions.cmdstan = str(c_beg+1:end-2);


%check the version of MatlabStan

%unfortunately, the release version of matlabstan is not contained in any
%of the data files. That means that I cannot check programmatically whether
%the version is up to date. I've posted an issue on github to request that
%the version be included in the files somewhere. 
%The toolbox automatically names the directory with the version number
%included, but checking this would only work if the user used the toolbox
%to download the dependents. 
versions.matlabstan = 'unknown'; %not useful

%check the version of MatlabProcessManager
%find the file that has the information about the version
mpm_filepath = which('processManager.m');

%load the file and pull the version information
mpm_file = fileread(mpm_filepath);
C = strsplit(mpm_file,'version = ');
str = C{2};
c_beg = strfind(str,'''');
versions.matlabprocessmanager = str(c_beg(1)+1:c_beg(2)-1);

%print information to command window
fprintf('\nMatlab version: %s', versions.matlab);
fprintf('\nERA Toolbox version: %s', versions.ERA);
fprintf('\nCmdStan version: %s', versions.cmdstan);
fprintf('\nMatlabStan version: %s', versions.matlabstan);
fprintf('\nMatlab Process Manager version: %s\n\n',...
    versions.matlabprocessmanager);

end

