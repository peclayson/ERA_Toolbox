function era_installdependents(varargin)
%Install the Matlab dependents. User will be prompted when input necessary
%
%era_installdependents
%
%Last Updated 9/11/20
%
%Required Inputs:
% No inputs are required.
% Some features only work for certain versions of Mac or Windows. When
%  necessary, the user will be given instructions re: installation of
%  various pieces by being directed to the User Manual.
%Optional Inputs:
% depvercheck - structure array with three dependents. The script will
%  update all of the dependents with a value of 0 (indicating that the
%  particular dependnet is out of date).
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
% by Peter Clayson (7/31/16)
% peter.clayson@gmail.com
%
%8/6/16 PC
% for some reason the untar doesn't properly unpack cmdstan on OS X 
%  Yosemite. Add check to install the .zip if using OS X Yosemite. 
% added gui for user to input when XCode installation is complete
% added more cw updates
% added capability to install dependents on Windows 7
%
%9/18/16 PC
% added paths to updated cmdstan
%
%11/10/16 PC
% added urls variable so I don't have to go through and manually update
%  each url when there is an update
%
%1/19/17 PC
% updated copyright
% updated cmdstan and matlab process manager versions
%
%6/22/17 PC
% updated dependents version
% minor changes
%
%6/24/17 PC
% version numbers are now pulled from the era_dependentsversions file
% software can now be updated with this script
%
%6/25/17 PC
% fix depvercheck when no inputs provided
% fixed some wording of the warnings
% made installation gui a bit wider
%
%8/15/17 PC
% fixed problem checking whether linux is being used
%
%9/8/17 PC
% Macs will not use the terminal command for unpacking the tarball. I was
%  getting a weird error when using Matlab's unpack for the new version of
%  CmdStan. Switching to the terminal command seems to have fixed it.
%
%6/17/19 PC
% Fixed bug where cmdstan path not updated in stan_home when cmdstan was 
%  updated on Mac OS
%
%9/11/20 PC
% Changes assocaited with fixing MatlabStan installation (need to manually
%  input version number into script
%
%1/20/21 PC
% Changes associated with checking the version of Mac OS >11
%
%3/5/21 PC
% fix problem with checking compiler on Windows
% fixed url to .zip for cmdstan
%
%3/18/21 PC
% continued fixes to Windows installation


%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if nargin > 0 && ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_installdependents for more information on optional inputs'));
    end
    
    %check if the dataset is present
    ind = find(strcmp('depvercheck',varargin),1);
    if ~isempty(ind)
        depvercheck = varargin{ind+1}; 
    else 
        depvercheck = [];
    end
elseif nargin ~= 2
    depvercheck = [];
end

%load the versions of the dependents used by the toolbox
depvers = era_dependentsversions;

url_cs_zip = strcat('https://github.com/stan-dev/cmdstan/archive/v',...
    depvers.cmdstan,'.zip');
url_cs_tar = strcat('https://github.com/stan-dev/cmdstan/releases/download/v',...
    depvers.cmdstan,'/cmdstan-',depvers.cmdstan,'.tar.gz');
url_ms = strcat('https://github.com/brian-lau/MatlabStan/archive/v',...
    depvers.matlabstan);
url_mp = strcat('https://github.com/brian-lau/MatlabProcessManager/archive/v',...
    depvers.matlabprocessmanager);

urls = struct;
urls.cmdstan_zip = url_cs_zip;
urls.cmdstan_tar = url_cs_tar;
urls.ms_zip = strcat(url_ms,'.zip');
urls.ms_tar = strcat(url_ms,'.tar.gz');
urls.mpm_zip = strcat(url_mp,'.zip');
urls.mpm_tar = strcat(url_mp,'.tar.gz');

%determine the version of OS that is being used
if ismac
    sys = 1; %because Apple is number 1, obviously
elseif ispc
    sys = 2; %because Windows OS is inferior
elseif isunix && ~ismac
    error('os:linux',... %Error code and associated error
    strcat('WARNING: Automatic installation does not work for Linux \n\n',...
    'See help era_installdependents for more information'));
end

%check whether the version of the OS is supported
if sys == 1 %mac
    [~,cmdout] = system('sw_vers');
    parseout = strsplit(cmdout);
    parseOS = strsplit(parseout{6},'.');
    
    if isnan(str2double(parseOS{1})) || ~contains(char(parseOS{1}),'.')
        parseOS = strsplit(parseout{4},'.');
    end
    
    if str2double(parseOS{1}) >= 10
        if str2double(parseOS{2}) >= 9 || str2double(parseOS{1}) >= 11
            
            if isempty(depvercheck)
            
                [~, savepath] = uiputfile('ERADependents',...
                'Where would you like to save the directory for the dependents?');

                %if the user does not select a file, then take the user back to era_start    
                if savepath == 0 
                    errordlg('Location not selected','File Error');
                    era_start;
                    return;
                end

                wrkdir = fullfile(savepath,'ERADependents');
                
            elseif ~isempty(depvercheck)
                
                wrkdir = fileparts(fileparts(which('runCmdStanTests.py')));
                
                rmolddeps(depvercheck);
                
            end
            
            era_macdepsinstall(wrkdir,str2double(parseOS{2}),urls);
            
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
    
    if isempty(depvercheck)
    
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
     
    elseif ~isempty(depvercheck)
                
        wrkdir = fileparts(fileparts(which('runCmdStanTests.py')));

        rmolddeps(depvercheck);
                
    end
    
    era_windepsinstall(wrkdir,urls);
end
        
end

function era_macdepsinstall(wrkdir,OSver,urls)

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
    loc = urls.cmdstan_tar;

    %change the working directory to where the files will be installed
    cd(wrkdir);

    %download the cmdstan tarball
    fileout = websave('cmdstan.tar.gz',loc);

    %unpack the tarball 
    %for some reason Matlab's untar function does not unpack the file
    %correctly, so the system command is used instead
    status = system('tar --no-same-owner -x -z -f cmdstan.tar.gz');    

    %verify that the tarball was successfully unpacked
    if status == 1
        error('CmdStan:notinstalled',... %Error code and associated error
            strcat('ERROR: The system failed to unpack the cmdstan tarball\n\n',...
            'Please verify that C++ development environment is installed\n',...
            'The recommended C++ development environment is XCode\n',...
            'For installation instructions, see Mac Installation Instructions\n',...
            'in Appendix A of the User Manual'));
    end

    %delete the tarball after it's been unpacked
    delete(fullfile(wrkdir,'cmdstan.tar.gz'));

    %now find the path to the new cmdstan directory
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);
    
    if length(ind>1) %#ok<*ISMT>
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    
    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan succesfully downloaded and unpacked\n');
    era_guistatus(depcheck);
    
elseif depcheck.cmdstan == 0 && OSver == 10
    %url for cmdstan tarball
    loc = urls.cmdstan_zip;

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
    
    if length(ind>1)
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    
    %update depcheck and display status
    depcheck.cmdstan = 1;
    fprintf('CmdStan succesfully downloaded and unzipped\n');
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
            strcat('ERROR: Automatic build of CmdStan failed \n\n',...
            'Please install CmdStan manually\n',...
            'For installation instructions, see Appendix A of the User Manual'));
    end
    
    newcmdstanbuild = 1;
    
end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %url for Matlab Process Manager
    loc = urls.mpm_tar;
    
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
    loc = urls.ms_tar;
    
    %download tarball
    fileout = websave('ms.tar.gz',loc);
    
    %unpack tarball
    untar(fileout,wrkdir);
    
    %get dir of MatlabStan so stan_home can be edited
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    
    if length(ind>1)
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
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
    
    
    %open the StanModel.m file and change the version to current cmdstan
    % version
    fid=fopen(fullfile(msdir,'StanModel.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    ind = strcmp(storefile,"            ver = cellfun(@str2num,regexp(str{3},'\.','split'));");
    
    jvers = era_dependentsversions;
    storefile{ind} = ['            ver = ''' jvers.cmdstan ''';'];

    %overwrite the existing stan_home.m file with the new information
    fid=fopen(fullfile(msdir,'StanModel.m'),'w');
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
  
%this is necessary in case cmdstan was updated, but a new version  
%matlabstan was not installed. Matlabstan needs to be pointed to cmdstan.
elseif depcheck.mstan == 1 && exist('newcmdstanbuild','var')
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);
    
    %get dir of MatlabStan so stan_home can be edited
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    
    if length(ind>1)
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    msdir = fullfile(wrkdir,ls(ind).name);

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


function era_windepsinstall(wrkdir,urls)
%install dependents on Windows Operating System

%get the current directory so it can be reverted to after installation
startdir = cd;

%create a structure array that will store which dependents are not properly
%installed
%0 - not installed, 1 - installed
depcheck = struct;

%first check whether g++ and make compilers exist
[~,cmdout] = system('make -v');
makematch = strncmp(cmdout,'GNU Make',8);

% makematch = 1;

if makematch == 1
    depcheck.CLT = 1;
    fprintf('Command line tools installed\n');
else
    error('Rtools:notinstalled',... %Error code and associated error
        strcat('WARNING: Command line tools not installed\n\n',...
        'Automatic installation of command line tools not supported for Windows\n',...
        'See User Manual for instructions on how to install the Rtools package\n',...
        'Rtools contains the necessary command line tools needed for CmdStan\n',...
        'Sorry I was unable to automate this process for Windows users!\n\n')); 
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

    %url for cmdstan zip
    loc = urls.cmdstan_zip;

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
    
    if length(ind>1)
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
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
    
    cmdstannewbuild = 1;
end

%check whether the Matlab Process Manager needs to be installed
if depcheck.mpm == 0
    %url for Matlab Process Manager
    loc = urls.mpm_zip;
    
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
    loc = urls.ms_zip;
    
    %download tarball
    fileout = websave('ms.tar.gz',loc);
    
    %unpack tarball
    unzip(fileout,wrkdir);
    
    %get dir of MatlabStan so stan_home can be edited
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    
    if length(ind>1)
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    msdir = fullfile(wrkdir,ls(ind).name);
    
    %delete the odd file that can be created some times
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end
    
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if ~exist('cmdstandir','var') || isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    elseif exist('cmdstandir','var') && isempty(cmdstandir)
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
    
    %open the StanModel.m file and change the version to current cmdstan
    % version
    fid=fopen(fullfile(msdir,'StanModel.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);

    ind = strcmp(storefile,"            ver = cellfun(@str2num,regexp(str{3},'\.','split'));");
    
    jvers = era_dependentsversions;
    storefile{ind} = ['            ver = ''' jvers.cmdstan ''';'];

    %overwrite the existing stan_home.m file with the new information
    fid=fopen(fullfile(msdir,'StanModel.m'),'w');
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
    
%this is necessary if cmdstan is updated. Matlabstan needs to have the 
%location of cmdstan updated.
elseif depcheck.mstan == 1 && exist('cmdstannewbuild','var')
    %in case cmdstan wasn't installed in this run, find the location of the
    %cmdstan dir
    if isempty(cmdstandir)
        p = which('runCmdStanTests.py');
        cmdstandir = fileparts(p);
    end
    cd(cmdstandir);
    
    %get dir of MatlabStan so stan_home can be edited
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    
    if length(ind>1)
        [~,newind] = max([ls(ind).datenum]);
        ind = ind(newind);
    end
    
    msdir = fullfile(wrkdir,ls(ind).name);

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
figwidth = 600;
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

function rmolddeps(depvercheck)

fprintf('Removing old versions of dependents\nThis may take a while\n\n');

%temporarily turn off warnings that filepaths are being removed
warning('off','MATLAB:rmpath:DirNotFound');
warning('off','MATLAB:RMDIR:RemovedFromPath');

if depvercheck.cmdstan == 0     
    %delete old cmdstan files
    rmdir(fileparts(which('runCmdStanTests.py')),'s'); 
end

if depvercheck.matlabstan == 0 
    %delete old matlabstan files
    rmdir(fileparts(which('mcmc.m')),'s');  
end

if depvercheck.matlabprocessmanager == 0 
    %delete old matlabprocessmanager files
    rmdir(fileparts(which('processManager.m')),'s');  
end

%turn warnings back on
warning('on','MATLAB:rmpath:DirNotFound');
warning('on','MATLAB:RMDIR:RemovedFromPath');

end






