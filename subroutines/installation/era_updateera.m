function era_updateera
%Update the ERA Toolbox. This currently only works for some newer versions
% of Mac or Windows. If using Linux, the ERA_Toolbox will need to be 
% updated manually.
%
%era_updateera
%
%Last Updated 6/17/17
%
%Required Inputs:
% No inputs are required.
% Some features only work for certain versions of Mac or Windows.
%
%Output:
% No data are outputted to the Matlab command window. However, software
%  will be installed and a new Matlab path will be saved.

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
% 10/6/17 PC
%  made changes associated with downloading the new .zip that is provided
%   so I can attempt to keep track of downloads
%
% 10/20/17 PC
%  accidentally left some code in that I used for debugging. It's been
%   taken out.

%link to newest ERA Toolbox release on github
urlraw = 'https://github.com/peclayson/ERA_Toolbox/releases/download/';
urlstr = 'https://github.com/peclayson/ERA_Toolbox/releases/latest';

webraw = webread(urlstr,'text','html');

%pull the string that contain the version number
rellist = regexp(webraw,'<h1 class="release-title">.*?</h1>','match');
str = strrep(rellist{1},'href','HREF');
str = strsplit(str,'>');

verraw = strsplit(str{strncmp('Version',str,7)},'<');
ver = strsplit(verraw{1},'Version ');

%version number will be appended to directory
ver = ver{2};

urlstr = strcat(urlraw,'v',ver);

era_dirname = strcat('ERA_Toolbox_v',ver);

str_cantdl = 'The new version of the toolbox is available on ';
str_cantdl = [str_cantdl... 
    '<a href="matlab:web(''https://github.com/peclayson/ERA_Toolbox/releases/latest'',''-browser'')">Github</a>'];

%determine the version of OS that is being used
if ismac
    sys = 1; %because Apple is number 1, obviously
elseif ispc
    sys = 2; %because Windows OS is inferior
elseif IsLinux
    error('os:linux',... %Error code and associated error
    strcat('WARNING: Automatic installation does not work for Linux \n\n',...
    'See help era_updateera for more information'));
end

%check whether the version of the OS is supported
if sys == 1 %mac
    [~,cmdout] = system('sw_vers');
    parseout = strsplit(cmdout);
    parseOS = strsplit(parseout{6},'.');
    if str2double(parseOS{1}) == 10
        if str2double(parseOS{2}) >= 9
            [~, savepath] = uiputfile('ERA_Toolbox',...
            'Where would you like to save the ERA_Toolbox?');
        
            %if the user does not select a file, then take the user back to era_start    
            if savepath == 0 
                errordlg('Location not selected','File Error');
                era_start;
                return;
            end
            
            wrkdir = fullfile(savepath,era_dirname);
            era_erainstall(wrkdir,urlstr);
            
        else
            error('mac:oldver',... %Error code and associated error
            strcat('WARNING: Automatic installation only works for OSX 10.9 and newer \n\n',...
            str_cantdl));
        end
    else
        error('mac:oldver',... %Error code and associated error
        strcat('WARNING: Automatic installation only works for Mac OS 10.9 or newer \n\n',...
            str_cantdl));
    end
elseif sys == 2 %windows
    fprintf('Warning: automatic installation has only been tested using Windows 7\n');
    [~, savepath] = uiputfile('ERA_Toolbox',...
    'Where would you like to save the directory for the updated ERA_Toolbox?');

    %if the user does not select a file, then take the user back to era_start    
    if savepath == 0 
        errordlg('Location not selected','File Error');
        era_start;
        return;
    end

    wrkdir = fullfile(savepath,era_dirname);

    era_erainstall(wrkdir,urlstr);
end
        
end

function era_erainstall(wrkdir,urlstr)
%install new version of ERA_Toolbox on a Mac

%get the current directory so it can be changed after installation
startdir = cd;

%get the directory of the current era_start file
old_eradir = which('era_start.m');
old_eradir = fileparts(old_eradir);

%add extenstion to html to download
urlstr = strcat(urlstr,'/ERA_Toolbox.zip');

%download zip
fileout = websave('era_toolbox.zip',urlstr);

%get the installation directory
installdir = fileparts(fileout);

%unzip
unzip(fileout,installdir);

%move the files to where they're supposed to go
movefile(fullfile(installdir,'ERA_Toolbox'),wrkdir)

%delete the zip file after it's been unpacked
delete(fileout);

%temporarily turn off warnings that filepaths are being removed
warning('off','MATLAB:rmpath:DirNotFound');

%remove old ERA_Toolbox files from the path
rmpath(genpath(old_eradir)); 

%turn warnings back on
warning('on','MATLAB:rmpath:DirNotFound');

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
savepath;
fprintf('Directories for ERA Toolbox added and saved to Matlab path\n');

%go back to the starting directory
if strcmp(genpath(old_eradir),startdir)
    cd(startdir);
else
    cd(wrkdir);
end

era_start;

end
