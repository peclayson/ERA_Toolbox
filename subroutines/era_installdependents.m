function era_installdependents
%Install the Matlab dependents. User will be prompted when input necessary
%
%era_installdependents
%
%Last Updated 8/6/16
%
%Required Inputs:
% No inputs are required.
% Some features only work for certain versions of Mac or Windows. When
%  necessary, the user will be given instructions re: installation of
%  various pieces by being directed to the User Manual.
%
%Output:
% No data are outputted to the Matlab command window. However, software
%  will be installed and a new Matlab path will be saved.

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
% by Peter Clayson (7/31/16)
% peter.clayson@gmail.com
%
%8/6/16 PC
% for some reason the untar doesn't properly unpack cmdstan on OS X 
%  Yosemite. Add check to install the .zip if using OS X Yosemite. 
% added gui for user to input when XCode installation is complete
% added more cw updates
% added capability to install dependents on Windows 7

%determine the version of OS that is being used
if ismac
    sys = 1; %because Apple is number 1, obviously
elseif ispc
    sys = 2; %because Windows OS is inferior
elseif IsLinux
    error('os:linux',... %Error code and associated error
    strcat('WARNING: Automatic installation does not work for Linux \n\n',...
    'See help era_installdependents for more information'));
end

%check whether the version of the OS is supported
if sys == 1 %mac
    [~,cmdout] = system('sw_vers');
    parseout = strsplit(cmdout);
    parseOS = strsplit(parseout{6},'.');
    if str2double(parseOS{1}) == 10
        if str2double(parseOS{2}) >= 9
            [~, savepath] = uiputfile('ERADependents',...
            'Where would you like to save the directory for the dependents?');
        
            %if the user does not select a file, then take the user back to era_start    
            if savepath == 0 
                errordlg('Location not selected','File Error');
                era_start;
                return;
            end
            
            wrkdir = fullfile(savepath,'ERADependents');
            era_macdepsinstall(wrkdir,str2double(parseOS{2}));
            
        else
            error('mac:oldver',... %Error code and associated error
            strcat('WARNING: Automatic installation only works for OSX 10.9 and newer \n\n',...
            'For installation instructions, see Appendix A of the User Manual'));
        end
    else
        error('mac:oldver',... %Error code and associated error
        strcat('WARNING: Automatic installation only works for Mac OS 10.9 or newer \n\n',...
        'For installation instructions, see Appendix A of the User Manual\n'));
    end
elseif sys == 2 %windows
    fprintf('Warning: automatic installation has only been tested using Windows 7\n');
    [~, savepath] = uiputfile('ERADependents',...
    'Where would you like to save the directory for the dependents?');

    %if the user does not select a file, then take the user back to era_start    
    if savepath == 0 
        errordlg('Location not selected','File Error');
        era_start;
        return;
    end

    wrkdir = fullfile(savepath,'ERADependents');
    
    era_windepsinstall(wrkdir);
end
        
end

function era_macdepsinstall(wrkdir,OSver)

%get the current directory so it can be reverted to after installation
startdir = cd;

%create a structure array that will store which dependents are not properly
%installed
%0 - not installed, 1 - installed
depcheck = struct;

%first check for the XCode command line tools
if exist('/Library/Developer/CommandLineTools','file') == 7
    depcheck.CLT = 1;
    fprintf('XCode Command line tools is installed\n');
else
    depcheck.CLT = 0;
    fprintf('XCode Command line tools does not appear to be installed\n');
end

%check for cmdstan
if exist('makefile','file') ~= 2 || ...
    exist('test-all.sh','file') ~= 2 || ...
    exist('runCmdStanTests.py','file') ~= 2 
   
    %CmdStan has not been installed
    depcheck.cmdstan = 0;
    
    %CmdStan has not been built
    depcheck.cmdbuild = 0;
    fprintf('CmdStan needs to be installed and built\n');
else
    
    %CmdStan has been installed
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');
    
    %check whether CmdStan has been properly built
    if exist('stansummary','file') ~= 2 || ...
        exist('stanc','file') ~= 2
   
        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end
end

%check for MatlabProcessManager (MPM)
if exist('processManager.m','file') ~= 2 || ...
    exist('processState.m','file') ~= 2 
    
    %MPM not installed
    depcheck.mpm = 0;
    fprintf('MatlabProcessManager is not installed\n');
else
    %MPM is installed
    depcheck.mpm = 1;
    fprintf('MatlabProcessManager is installed\n');
end

%check for MatlabStan
if exist('StanFit.m','file') ~= 2 || ...
    exist('StanModel.m','file') ~= 2 
    
    %mstan not installed
    depcheck.mstan = 0;
    fprintf('MatlabStan is not installed\n');
else
    %mstan is installed
    depcheck.mstan = 1;
    fprintf('MatlabStan is installed\n');
end

%pullup a gui to display the stauts of the ERA toolbox dependents
era_guistatus(depcheck);

%check whether the directory to save the files exists, if not create it
if exist(wrkdir,'file') ~=7
    %make the directory where the files will be saved
    mkdir(wrkdir);
end

%if needed, install CLT
if depcheck.CLT == 0
    %user will need to interacte with gui for installation of CLT
    status = system('xcode-select --install');
    
    fsize = get(0,'DefaultTextFontSize') + 3;
    figwidth = 400; figheight = 200;
    %initialize gui
    f = figure('unit','pix','Visible','off',...
      'position',[400 400 figwidth figheight],...
      'menub','no',...
      'numbertitle','off',...
      'resize','off');

    movegui(f,'center');

    %Write text
    uicontrol(f,'Style','text','fontsize',fsize+2,...
        'HorizontalAlignment','center',...
        'String','Please indicate when the XCode installation is complete',...
        'Position',[0 100 figwidth figheight/3]);          

    %Create a button that will resume installations
    uicontrol(f,'Style','push','fontsize',fsize,...
        'HorizontalAlignment','center',...
        'String','<html><center>Installation <br>Complete',...
        'Position', [(figwidth/5)/2 25 (figwidth/5)*4 75],...
        'Callback','uiresume(gcbf)'); 
    
    %display gui
    set(f,'Visible','on');
    %tag gui
    f.Tag = 'f';
    
    uiwait(f);
    
    %check if era_gui is open.
    f = findobj('Tag','f');
    if ~isempty(f)
        close(f);
    end
    
    if status == 0 
        depcheck.CLT = 1;
        fprintf('XCode Command Line Tools succesfully installed\n');
        %pullup a gui to display the stauts
        era_guistatus(depcheck);
    else
        error('CLT:notinstalled',... %Error code and associated error
        strcat('WARNING: Automatic installation of XCode CLT failed \n\n',...
        'Please install the command line tools manually\n',...
        'For installation instructions, see Appendix A of the User Manual'));
    end
end

%check whether cmdstan needs to be installed
if depcheck.cmdstan == 0 && OSver ~= 10

    %url for cmdstan tarball
    loc = 'https://github.com/stan-dev/cmdstan/releases/download/v2.11.0/cmdstan-2.11.0.tar.gz';

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.tar.gz',loc);

    %unpack the tarball
    untar(fileout,wrkdir);

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.tar.gz'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    
    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan succesfully installed\n');
    era_guistatus(depcheck);
    
elseif depcheck.cmdstan == 0 && OSver == 10
    %url for cmdstan tarball
    loc = 'https://github.com/stan-dev/cmdstan/releases/download/v2.11.0/cmdstan-2.11.0.zip';

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.zip',loc);

    %unpack the tarball
    unzip(fileout,wrkdir);

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.zip'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    
    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan succesfully installed\n');
    era_guistatus(depcheck);
end

%check whether CmdStan needs to be built
if depcheck.cmdbuild == 0
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var')
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    
    %change the cd to where cmdstan is installed so it can be built
    cd(cmdstandir);
    
    %find out how many cores are available for building stan
    ncores = feature('numCores');
    if ncores > 2
        ncores = ncores - 1;
    end

    %execute command to build stan
    buildcmdstan = strcat('make build -j',num2str(ncores));
    status = system(buildcmdstan);

    %update depcheck and display status
    if status == 0
        depcheck.cmdbuild = 1;
        fprintf('CmdStan successfully built\n');
        era_guistatus(depcheck);
    else
        error('CmdStan:notbuilt',... %Error code and associated error
            strcat('WARNING: Automatic build of CmdStan failed \n\n',...
            'Please install CmdStan manually\n',...
            'For installation instructions, see Appendix A of the User Manual'));
    end
end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %url for Matlab Process Manager
    loc = 'https://github.com/brian-lau/MatlabProcessManager/archive/v0.4.2.tar.gz';
    
    %download tarball
    fileout = websave('mpm.tar.gz',loc);
    
    %unpack tarball
    untar(fileout,wrkdir);
    
    %change status of depcheck and display gui
    depcheck.mpm = 1;
    fprintf('Matlab Process manager successfully installed\n');
    era_guistatus(depcheck);
end

%check whether MatlabStan is installed
if depcheck.mstan == 0
    %url for MatlabStan
    loc = 'https://github.com/brian-lau/MatlabStan/archive/v2.7.0.0.tar.gz';
    
    %download tarball
    fileout = websave('ms.tar.gz',loc);
    
    %unpack tarball
    untar(fileout,wrkdir);
    
    %get dir of MatlabStan so stan_home can be edited
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    msdir = fullfile(wrkdir,ls(ind).name);
    
    %delete the odd file that can be created some times
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end
    
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);

    %open the stan_home.m file and change the path to cmdstan
    fid=fopen(fullfile(msdir,'+mstan','stan_home.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    storefile{8} = ['d = ''' cmdstandir ''';'];

    %overwrite the existing stan_home.m file with the new information
    fid = fopen(fullfile(msdir,'+mstan','stan_home.m'),'w');
    for i=1:(length(storefile)-1)
        storefile{i} = strrep(storefile{i},'%','%%');
        storefile{i} = strrep(storefile{i},'\','\\');
        fprintf(fid,[storefile{i} '\n']);
    end
    fclose(fid);
    
    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan succesfully installed\n');
    era_guistatus(depcheck);
end

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
savepath;
fprintf('Directories for dependents added to Matlab path\n');

%close gui
era_instgui = findobj('Tag','era_instgui');

if ~isempty(era_instgui)
    
    close(era_instgui);
    
end

%go back to the starting directory
cd(startdir);

%rerun ERA Toolbox
era_start;

end


function era_windepsinstall(wrkdir)
%install dependents on Windows Operating System

%get the current directory so it can be reverted to after installation
startdir = cd;

%create a structure array that will store which dependents are not properly
%installed
%0 - not installed, 1 - installed
depcheck = struct;

%first check whether g++ and make compilers exist
[~,cmdout] = system('make -v');
makematch = strncmp(cmdout,'GNU Make version',16);

if makematch == 1
    depcheck.CLT = 1;
    fprintf('Command line tools installed\n');
else
    error('Rtools:notinstaled',... %Error code and associated error
        strcat('WARNING: Command line tools not installed\n\n',...
        'Automatic installation of command line tools not supported for Windows\n',...
        'See User Manual for instructions on how the Rtools package\n',...
        'Rtools contains the necessary command line tools needed for CmdStan\n',...
        'Sorry\n\n')); 
end

%check for cmdstan
if exist('makefile','file') ~= 2 || ...
    exist('test-all.sh','file') ~= 2 || ...
    exist('runCmdStanTests.py','file') ~= 2 
   
    %CmdStan has not been installed
    depcheck.cmdstan = 0;
    
    %CmdStan has not been built
    depcheck.cmdbuild = 0;
    fprintf('CmdStan needs to be installed and built\n');
else
    
    %CmdStan has been installed
    depcheck.cmdstan = 1;
    fprintf('CmdStan is installed\n');
    
    %check whether CmdStan has been properly built
    if exist('stansummary.exe','file') ~= 2 || ...
        exist('stanc.exe','file') ~= 2
   
        %CmdStan has not been built
        depcheck.cmdbuild = 0;
        fprintf('CmdStan has not been properly built\n');
    else
        %CmdStan has been built
        depcheck.cmdbuild = 1;
        fprintf('CmdStan has been properly built\n');
    end
end

%check for MatlabProcessManager (MPM)
if exist('processManager.m','file') ~= 2 || ...
    exist('processState.m','file') ~= 2 
    
    %MPM not installed
    depcheck.mpm = 0;
    fprintf('MatlabProcessManager is not installed\n');
else
    %MPM is installed
    depcheck.mpm = 1;
    fprintf('MatlabProcessManager is installed\n');
end

%check for MatlabStan
if exist('StanFit.m','file') ~= 2 || ...
    exist('StanModel.m','file') ~= 2 
    
    %mstan not installed
    depcheck.mstan = 0;
    fprintf('MatlabStan is not installed\n');
else
    %mstan is installed
    depcheck.mstan = 1;
    fprintf('MatlabStan is installed\n');
end

%pullup a gui to display the stauts of the ERA toolbox dependents
era_guistatus(depcheck);

%check whether the directory to save the files exists, if not create it
if exist(wrkdir,'file') ~=7
    %make the directory where the files will be saved
    mkdir(wrkdir);
end

%check whether cmdstan needs to be installed
if depcheck.cmdstan == 0 

    %url for cmdstan tarball
    loc = 'https://github.com/stan-dev/cmdstan/releases/download/v2.11.0/cmdstan-2.11.0.zip';

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.zip',loc);

    %unpack the tarball
    unzip(fileout,wrkdir);

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.zip'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    
    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan succesfully installed\n');
    era_guistatus(depcheck);
    
end

%check whether CmdStan needs to be built
if depcheck.cmdbuild == 0
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var')
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    
    %change the cd to where cmdstan is installed so it can be built
    cd(cmdstandir);
    
    %find out how many cores are available for building stan
    ncores = feature('numCores');
    if ncores > 2
        ncores = ncores - 1;
    end

    %execute command to build stan
    buildcmdstan = strcat('make build -j',num2str(ncores));
    status = system(buildcmdstan);

    %update depcheck and display status
    if status == 0
        depcheck.cmdbuild = 1;
        fprintf('CmdStan successfully built\n');
        era_guistatus(depcheck);
    else
        error('CmdStan:notbuilt',... %Error code and associated error
            strcat('WARNING: Automatic build of CmdStan failed \n\n',...
            'Please install CmdStan manually\n',...
            'For installation instructions, see Appendix A of the User Manual'));
    end
end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %url for Matlab Process Manager
    loc = 'https://github.com/brian-lau/MatlabProcessManager/archive/v0.4.2.zip';
    
    %download tarball
    fileout = websave('mpm.zip',loc);
    
    %unpack tarball
    unzip(fileout,wrkdir);
    
    %change status of depcheck and display gui
    depcheck.mpm = 1;
    fprintf('Matlab Process manager successfully installed\n');
    era_guistatus(depcheck);
end

%check whether MatlabStan is installed
if depcheck.mstan == 0
    %url for MatlabStan
    loc = 'https://github.com/brian-lau/MatlabStan/archive/v2.7.0.0.zip';
    
    %download tarball
    fileout = websave('ms.tar.gz',loc);
    
    %unpack tarball
    unzip(fileout,wrkdir);
    
    %get dir of MatlabStan so stan_home can be edited
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    msdir = fullfile(wrkdir,ls(ind).name);
    
    %delete the odd file that can be created some times
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end
    
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);

    %open the stan_home.m file and change the path to cmdstan
    fid=fopen(fullfile(msdir,'+mstan','stan_home.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    storefile{8} = ['d = ''' cmdstandir ''';'];

    %overwrite the existing stan_home.m file with the new information
    fid = fopen(fullfile(msdir,'+mstan','stan_home.m'),'w');
    for i=1:(length(storefile)-1)
        storefile{i} = strrep(storefile{i},'%','%%');
        storefile{i} = strrep(storefile{i},'\','\\');
        fprintf(fid,[storefile{i} '\n']);
    end
    fclose(fid);
    
    %update depcheck and display gui
    depcheck.mstan = 1;
    fprintf('MatlabStan succesfully installed\n');
    era_guistatus(depcheck);
end

%add files to the Matlab path and save it
addpath(genpath(wrkdir));
savepath;
fprintf('Directories for dependents added to Matlab path\n');

%close gui
era_instgui = findobj('Tag','era_instgui');

if ~isempty(era_instgui)
    
    close(era_instgui);
    
end

%go back to the starting directory
cd(startdir);

%rerun ERA Toolbox
era_start;

end


function era_guistatus(depcheck)

%define parameters for figure position
figwidth = 400;
figheight = 450;
fsize = get(0,'DefaultTextFontSize') + 3;

%define space between rows and first row location
rowspace = 35;
row = figheight - rowspace*2;

%define locations of column 1 and 2
lcol = 30;
rcol = (figwidth/8)*5;

era_instgui = findobj('Tag','era_instgui');

if isempty(era_instgui)

    %create the gui
    era_instgui= figure('unit','pix','Visible','off',...
      'position',[400 400 figwidth figheight],...
      'menub','no',...
      'name','Installation Summary',...
      'numbertitle','off',...
      'resize','off');

    movegui(era_instgui,'center');
    
end

%Print the name of the gui
uicontrol(era_instgui,'Style','text','fontsize',fsize+4,...
    'HorizontalAlignment','center',...
    'String','Installation Status',...
    'Position',[0 row figwidth 25]);     

row = row - rowspace*2;

%CLI status
uicontrol(era_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Command Line Tools',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.CLT == 0
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.CLT == 1
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%CmdStan status
uicontrol(era_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','CmdStan Files',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.cmdstan == 0
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.cmdstan == 1
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%CmdStan build status
uicontrol(era_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','CmdStan Build',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.cmdbuild == 0
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Not Built',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.cmdbuild == 1
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Built',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%Matlab Process Manager status
uicontrol(era_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Matlab Process Manager',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.mpm == 0
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.mpm == 1
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace;

%mstan status
uicontrol(era_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','left',...
    'String','Matlab Stan',...
    'Position', [lcol row figwidth/3 25]);  

if depcheck.mstan == 0
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Not Installed',...
        'Position', [rcol row figwidth/4 25]);  
elseif depcheck.mstan == 1
    uicontrol(era_instgui,'Style','text','fontsize',fsize,...
        'String','Installed',...
        'Position', [rcol row figwidth/4 25]);  
end

row = row - rowspace*2;

%mstan status
uicontrol(era_instgui,'Style','text','fontsize',fsize,...
    'HorizontalAlignment','center',...
    'String','Installation may take a while...',...
    'Position',[0 row figwidth 25]);

%make sure the gui is displayed if it has not been shown. Tag the gui as
%well
if strcmp(get(era_instgui,'Visible'),'off')
    set(era_instgui,'Visible','on');
    
    %tag gui
    era_instgui.Tag = 'era_instgui';
end

%pause for the gui to be displayed
pause(.02);

end









