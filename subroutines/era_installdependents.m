function era_installdependents

%test script to install cmdstan

%!xcodebuild -version
%

%check whether command line tools is installed. status = 0, if installed,
%else means not installed
%[status,cmdout] = system('xcode-select -p');
%[status,~] = system('xcode-select -p');

status = exist('/Library/Developer/CommandLineTools','file');

if status == 7
    fprintf('XCode Command line tools is installed\n');
else
    fprintf('XCode Command line tools does not appear to be installed\n');
end

%cmdstan
if exist('makefile','file') ~= 2 || ...
        exist('test-all.sh','file') ~= 2 || ...
        exist('runCmdStanTests.py','file') ~= 2 
else
    loc = 'https://github.com/stan-dev/cmdstan/releases/download/v2.11.0/cmdstan-2.11.0.tar.gz';
    wrkdir = '/Users/Peter/Documents/MATLAB/testdownload/';
    cd(wrkdir);
    fileout = websave('cmdstan.tar.gz',loc);
    untar(fileout,wrkdir);
    delete([wrkdir 'cmdstan.tar.gz']);
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'cmdstan', 7)==1);
    cmdstandir = fullfile(wrkdir,ls(ind).name);
    cd(cmdstandir);
    
    ncores = feature('numCores');
    if ncores > 2
        ncores = ncores - 1;
    end
    
    buildcmdstan = strcat('make build -j',num2str(ncores));
    [status,cmdout] = system(buildcmdstan);
    
    loc = 'https://github.com/brian-lau/MatlabProcessManager/archive/v0.4.2.tar.gz';
    fileout = websave('mpm.tar.gz',loc);
    untar(fileout,wrkdir);
    %delete([wrkdir 'mpm.tar.gz']);
    
    loc = 'https://github.com/brian-lau/MatlabStan/archive/v2.7.0.0.tar.gz';
    fileout = websave('ms.tar.gz',loc);
    untar(fileout,wrkdir);
    %delete([wrkdir 'ms.tar.gz']);
    ls = dir(wrkdir);
    ind = find(strncmp({ls.name}, 'MatlabStan', 10)==1);
    msdir = fullfile(wrkdir,ls(ind).name);
    cd(cmdstandir);
    
    if exist(fullfile(wrkdir,'pax_global_header'),'file') ~= 0
        delete(fullfile(wrkdir,'pax_global_header'));
    end
    
    fid=fopen(fullfile(msdir,'+mstan','stan_home.m'));
    storefile = {};
    while 1
        tline = fgetl(fid);
        storefile{end+1} = tline;
        if ~ischar(tline), break, end
    end
    fclose(fid);
        
    storefile{8} = ['d = ''' cmdstandir ''';'];

    fid = fopen(fullfile(msdir,'+mstan','stan_home.m'),'w');
    for i=1:(length(storefile)-1)
        storefile{i} = strrep(storefile{i},'%','%%');
        storefile{i} = strrep(storefile{i},'\','\\');
        fprintf(fid,[storefile{i} '\n']);
    end
    fclose(fid);

    addpath(genpath(wrkdir));
    savepath;
end







end


function era_guistatus(str)










end









