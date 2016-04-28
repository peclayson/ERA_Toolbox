function era_updatecheck(eraver)
%
%Check whether a new release has been posted on Github
%
%era_updatecheck
%
%Lasted Updated 4/27/16
%
%Required Input:
% eraver - ERA Toolbox version
%
%Outputs:
% No outputs. 
% Information will be printed to the command window for the user
%

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
% by Peter Clayson (4/27/16)
% peter.clayson@gmail.com

try
    %pull webpage from github
    urlstr = 'https://github.com/peclayson/ERA_Toolbox/releases';
    webraw = webread(urlstr,'text','html');

    %pull the string that contain the version number
    rellist = regexp(webraw,'<h1 class="release-title">.*?</h1>','match');
    str = strrep(rellist{1},'href','HREF');
    str = strsplit(str,'>');

    verraw = strsplit(str{strncmp('Version',str,7)},'<');
    ver = strsplit(verraw{1},'Version ');
    ver = ver{2};
    
    %let user know what was found
    if strcmp(ver,eraver)
        str = 'You are running the most up-to-date version of the toolbox';
        fprintf('\n%s\n',str);
    else
        str = 'There is a new version of the toolbox available on ';
        str = [str... 
            '<a href="matlab:web(''https://github.com/peclayson/ERA_Toolbox/releases'',''-browser'')">Github</a>'];
        fprintf('\n%s\n',str);
    end
catch
    fprintf('\nUnable to connect to Github to check for new releases\n');
end


end
