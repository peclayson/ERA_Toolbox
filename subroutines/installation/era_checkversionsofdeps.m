function versions = era_checkversionsofdeps
%Check the versions of the dependents for ERA Toolbox
%
%Last Updated 9/10/20
%

%This function just will examine whether the dependents are up to date. It
%currently does not work for MatlabStan because the version info is not
%contained in any of the downloaded files.
%
%Input
% No inputs required in the command line
%
%Output
% versions - array with information regarding whether versions of the
% dependents are 0-old, 1-up to date, or 2-untested new versions
%   cmdstan - information for cmdstan
%   matlabstan - information for matlab stan
%   matlabprocessmanager - information for matlab process manager
%

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
% by Peter Clayson (6/24/17)
% peter.clayson@gmail.com
%
%9/10/20 Pc
% changes for checking version of cmdstan after updating to 2.24.1

%create an empty array for storing the information about dependents
versions = struct;
versions.cmdstan = [];
versions.matlabstan = [];
versions.matlabprocessmanager = [];

%pull the current versions of the dependents used by the toolbox
depvers = era_dependentsversions;

%start off checking the version of cmdstan
%find the file that has the information about the version
%need to do this two different ways due to changes around ~vcmdstan 2.22
cs_filepath = which('cmdstan-guide.tex');

if ~isempty(cs_filepath)
    %load the file and pull the version information
    cs_file = fileread(cs_filepath);
    C = strsplit(cs_file,'\');
    str = C{5};
    c_beg = strfind(str,'{');
    cs_ver = str(c_beg+1:end-2);
    
else
    cs_filepath = fileparts(which('runCmdStanTests.py'));
    command = [fullfile(cs_filepath,'bin','stanc') ' --version'];
    %command = [fullfile(self.stan_home,'bin','stanc') ' --version'];
    p = processManager('id','stanc version','command',command,...
        'keepStdout',true,...
        'printStdout',false,...
        'pollInterval',0.005);
    
    p.block(0.05);
    
    str = regexp(p.stdout{1},'\ ','split');
    cs_ver = str{2}(2:end);
end

%format version information for cmdstan that is used by the toolbox
cs_used_parts = sscanf(depvers.cmdstan,'%d.%d.%d')';

%format version information for cmdstan that is found in the matlab path
cs_found_parts = sscanf(cs_ver,'%d.%d.%d')';

%compare cmdstan versions
for ii = 1:3
    if cs_used_parts(ii) > cs_found_parts(ii)
        versions.cmdstan = 0;
        break;
    elseif cs_used_parts(ii) < cs_found_parts(ii)
        versions.cmdstan = 2;
        break;
    elseif cs_used_parts(ii) == cs_found_parts(ii)
        versions.cmdstan = 1;
    end
end


%check the version of MatlabStan

%unfortunately, the release version of matlabstan is not contained in any
%of the data files. That means that I cannot check programmatically whether
%the version is up to date. I've posted an issue on github to request that
%the version be included in the files somewhere. 
%The toolbox automatically names the directory with the version number
%included, but checking this would only work if the user used the toolbox
%to download the dependents. 
versions.matlabstan = -1; %not useful


%check the version of MatlabProcessManager
%find the file that has the information about the version
mpm_filepath = which('processManager.m');

%load the file and pull the version information
mpm_file = fileread(mpm_filepath);
C = strsplit(mpm_file,'version = ');
str = C{2};
c_beg = strfind(str,'''');
mpm_ver = str(c_beg(1)+1:c_beg(2)-1);

%format version information for cmdstan that is used by the toolbox
mpm_used_parts = sscanf(depvers.matlabprocessmanager,'%d.%d.%d')';

%format version information for cmdstan that is found in the matlab path
mpm_found_parts = sscanf(mpm_ver,'%d.%d.%d')';

%compare cmdstan versions
for ii = 1:3
    if mpm_used_parts(ii) > mpm_found_parts(ii)
        versions.matlabprocessmanager = 0;
        break;
    elseif mpm_used_parts(ii) < mpm_found_parts(ii)
        versions.matlabprocessmanager = 2;
        break;
    elseif mpm_used_parts(ii) == mpm_found_parts(ii)
        versions.matlabprocessmanager = 1;
    end
end


end

