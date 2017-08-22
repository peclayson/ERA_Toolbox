function era_data = era_computerelwrap(varargin)
%
%Function to parse inputs from  into era_computerel
%
%era_dataout = era_relwrap('era_prefs',era_prefs,'era_data',era_data)
%
%Last Updated 8/22/17
%
%Input
% era_prefs - toolbox preferences
% era_data - toolbox data
%
%Output
% era_data - ERA Toolbox dataset with the variance components and
%  convergence checks from Stan in era_data.rel
% A directory will be created where the temporary files will be saved for
%  running the Stan model (stan creates various temp files). After the
%  model is done running (executed with era_computerel), the temporary
%  directory and files will be removed.

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
% by Peter Clayson (7/22/16)
% peter.clayson@gmail.com
%
%7/26/16 PC
% bug fix: not able to execute era_startproc_gui directly
%
%7/28/16 PC
% fixed typo causing era_start not to execute after failing to converge
%
%1/19/17 PC
% updated copyright
%
%8/15/17 PC
% fixed typos in comments
%
%8/16/17 PC
% added code to view trace plots prior to saving Stan outputs
%
%8/22/17 PC
% updated bug fixes to viewing trace plots

%pull era_prefs and era_data from varargin
[era_prefs, era_data] = era_findprefsdata(varargin);

%In case era_computerelwrap was executed using the CLI, check to make sure
%the data are set up.

%check for era_prefs
if isempty(era_prefs)
    error('era_prefs:notfound',...
        strcat('ERROR: Inputs are incomplete \n\n',...
        'era_prefs not found.\n',...
        'See help for era_computerelwrap more information'));
end

%check for era_data
if isempty(era_data)
    error('era_data:notfound',...
        strcat('ERROR: Inputs are incomplete \n\n',...
        'era_data not found.\n',...
        'See help for era_computerelwrap more information'));
end
  
%check that data are loaded
if ~isfield(era_data.proc,'data')
    %if the data have not been loaded where they are expected, try loading
    %the file specified in era_data
    try
        era_data.proc.data = era_loadfile('era_prefs',era_prefs,...
        'era_data',era_data);
    catch
        error('era_data:datanotfound',...
            strcat('ERROR: Inputs are incomplete\n\n',...
            'loaded dataset not found in era_data.proc.data\n',...
            'See help for era_loadfile to load dataset'));
    end
end

%check if there is a name of the dataset to be saved
if ~isfield(era_data.proc,'savename')
    error('era_data:savename',...
        strcat('ERROR: Inputs are incomplete\n\n',...
        'the name of the file to be saved has not been specified\n',...
        'information should be in era_data.proc.savename'));
else
    %ensure something is in savename
    if isempty(era_data.proc.savename)
        error('era_data:savename',...
            strcat('ERROR: Inputs are incomplete\n\n',...
            'the name of the file to be saved has not been specified\n',...
            'information should be in era_data.proc.savename'));
    end
    %check for whitespace for cmdstan
    if any(isspace(era_data.proc.savename))
        error('era_data:savename',...
            strcat('ERROR: Inputs are incomplete\n\n',...
            'cmdstan does not accept filenames with whitespace\n',...
            'please specify a filename without whitespace'));
    end
end

%check if the path for the new dataset exists
if ~isfield(era_data.proc,'savepath')
    error('era_data:savepath',...
        strcat('ERROR: Inputs are incomplete\n\n',...
        'the save path has not been specified\n',...
        'information should be in era_data.proc.savepath'));
else
    %check if the path is in savepath
    if isempty(era_data.proc.savepath)
        error('era_data:savepath',...
            strcat('ERROR: Inputs are incomplete\n\n',...
            'the directory for savepath does not exist\n',...
            'specified path:',era_data.proc.savepath));
    end
    %check if the path exists
    if ~exist(era_data.proc.savepath,'dir')
        error('era_data:savepath',...
            strcat('ERROR: Inputs are incomplete\n\n',...
            'the directory for savepath does not exist\n',...
            'specified path:',era_data.proc.savepath));
    end
    %check for whitespace for cmdstan
    if any(isspace(era_data.proc.savename))
        error('era_data:savepath',...
            strcat('ERROR: Inputs are incomplete\n\n',...
            'cmdstan does not accept paths with whitespace\n',...
            'please specify a path without whitespace'));
    end
end

%Change working dir for temporary Stan files
if ~exist(fullfile(era_data.proc.savepath,'Temp_StanFiles'), 'dir')
  mkdir(era_data.proc.savepath,'Temp_StanFiles');
end
origdir = cd(fullfile(era_data.proc.savepath,'Temp_StanFiles'));

%set initial state of rerun to 1
%this will run era_computerel
%if chains properly converged then rerun will be changed to 0 and the while
%loop will be exited
rerun = 1;

%whether chains converged will be checked each time era_computerel is run
while rerun ~= 0

    %pass the data to era_computerel for analysis
    REL = era_computerel('data',era_data.proc.data,...
        'chains',era_prefs.proc.nchains,...
        'iter',era_prefs.proc.niter,...
        'verbose',era_prefs.proc.verbose,...
        'showgui',2);
    
    %check convergence of chains
    RELout = era_checkconv(REL);
    
    if RELout.out.conv.converged == 0
       
        era_reruncheck;
        era_gui = findobj('Tag','era_gui');
        rerun = guidata(era_gui);
        close(era_gui);
        delete(era_gui);
        
        era_tplots = findobj('Tag','tplots');
        if ~isempty(era_tplots)
            close(era_tplots);
            delete(era_tplots);
        end
    else
        
        %if chains converged, do not rerun
        rerun = 0;
        fprintf('\nModel converged\n');
        
    end
    
    %check if the user wanted to see the trace plots
    if era_prefs.proc.traceplots == 2 && rerun == 0
        era_checktraceplots(REL,'askuser',1);
        era_gui = findobj('Tag','era_gui');
        rerun = guidata(era_gui);
        close(era_gui);
        delete(era_gui);
    end
    
    %if convergence was not met and the user would like to rerun the model,
    %double the number of iterations (if the doubled number is less than
    %1000, then run 1000 iterations)
    if rerun == 1
        era_prefs.proc.niter = era_prefs.proc.niter * 2;
        if era_prefs.proc.niter < 1000
            era_prefs.proc.niter = 1000;
        end
        fprintf('\nIncreasing number of iterations to %d\n',...
            era_prefs.proc.niter);
    end
end

%sometimes the era_gui doesn't close
era_gui = findobj('Tag','era_gui');
if ~isempty(era_gui)
    close(era_gui);
end
    
%change the working directory back to the original directory 
cd(origdir);

%attempt to remove the temporary directory and its files
try
    
    rmdir(fullfile(era_data.proc.savepath,'Temp_StanFiles'),'s');

catch
    
    fprintf('\n\nTemporary directory could not be removed.');
    fprintf('\nPath: %s\n',era_data.proc.savepath);
end

%if the chains did not converge send the user back to era_startproc_gui
if RELout.out.conv.converged == 0

    era_startproc('era_prefs',era_prefs,'era_data',era_data);
    return;
    
else
    
    %chains converged. Let the user know.
    %str = 'Models successfully converged with %d chains and %d iterations';
    %fprintf(strcat('\n',str,'\n'),procpref.nchains,procprefs.niter);
    
end

fprintf('\nSaving Processed Data...\n\n');

era_data.rel = RELout;

save(fullfile(era_data.proc.savepath,era_data.proc.savename),'era_data');

end
