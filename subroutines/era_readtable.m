function dataout = era_readtable(varargin)
%Load data into table by cycling through possible delimiters
%
%era_readtable('file','J:\Data\MyData.xlsx')
%
%Last Updated 1/19/17
%
%
%Input
% file - path and name of file to be loaded 
%  supported file extensions include .txt, .dat, .csv, .xls, .xlsx, .xlsb, 
%   .xlsm, .xltm, .xltx, and .ods
%
%Optional Input
% delimiter - delimiter used in text file to separate columns
%
%Output
% dataout - data in table format
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
% by Peter Clayson (4/18/16)
% peter.clayson@gmail.com
%
%1/19/17 PC
% updated copyright

%somersault through varargin inputs to check for which inputs were
%defined and store those values. 
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_readtable for more information on inputs'));
    end
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('file',varargin),1);
    if ~isempty(ind)
        file = cell2mat(varargin(ind+1)); 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end

    %check if a delimiter was specified 
    ind = find(strcmp('idcol',varargin),1);
    if ~isempty(ind)
        delimiter = cell2mat(varargin(ind+1)); 
    else
        delimiter = [];
    end
    
elseif ~isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help era_readtable for more information on inputs'));
    
end %if ~isempty(varargin)

%pull file extension 
[~,~,ext] = fileparts(file);

%file extensions able to be read by Matlab's readtable function
supportedfiles = {'txt','dat','csv','xls','xlsx','xlsb','xlsm','xltm',...
    'xltx','ods'};

%check if file is a supported file type
if ~isempty(ext) && strcmp(ext(1),'.')
    ext = ext(2:end);
end

%output error if file type is not supported
if sum(strcmp(ext,supportedfiles)) ~= 1

    error('ext:filetype',... %Error code and associated error
    strcat('WARNING: File type not supported \n\n',... 
    'For a list of supported file types, see help era_readtable \n'));

end

%turn off warning regarding Matlab reformatting headers
warning('off','MATLAB:table:ModifiedVarnames');

%set up a variable to specify whether the file has been loaded successfully
loadsuccess = 0;

%see if Matlab's readtable function can load the datafile simply based on
%the extension
if isempty(delimiter)
    dataout = readtable(file);
else
    dataout = readtable(file,'delimiter',delimiter);
end

%assume that if the data table has more than one column then it
%successfully loaded the data
if width(dataout) > 1
    loadsuccess = 1;
end

%make a variable containing some of the possible, common delimiters that 
%readtable can handle
possdel = {',',' ', '\t',';','|'};

%cycle through possible delimiters if the data file was not loaded
%successfully
i=1;
while loadsuccess == 0
    try 
        dataout = readtable(file,'delimiter',possdel{i});
        if width(dataout) > 1
            loadsuccess = 1;
        else
            if i == length(possdel) %stop cycling through while loop if  
                %there are no more delimiters to check
                loadsuccess = 3;
            end
            i = i+1;
        end
    catch
        if i == length(possdel) %stop cycling through while loop if there 
            %are no more delimiters to check
            loadsuccess = 3;
        end
        i = i+1;
    end
end

%turn (back) on warning regarding Matlab reformatting headers
warning('on','MATLAB:table:ModifiedVarnames');

if loadsuccess == 3
    error('loadingfile:delimiter',... %Error code and associated error
    strcat('WARNING: File type not successfully loaded \n\n',... 
    'For a list of supported file types, see help era_readtable \n'));
end

end

