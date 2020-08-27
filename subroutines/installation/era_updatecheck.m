function erainstalledver = era_updatecheck(eraver)
%Check whether a new release has been posted on Github
%
%era_updatecheck
%
%Lasted Updated 8/5/19
%
%Required Input:
% eraver - ERA Toolbox version
%
%Outputs:
% erainstalledver - contains information regarding whether the toolbox is
%  up to date 0 - old, 1 - current, 2 - beta
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
% by Peter Clayson (4/27/16)
% peter.clayson@gmail.com
%
%7/19/16 PC
% chaanged url to check 
%  https://github.com/peclayson/ERA_Toolbox/releases/latest
%
%1/19/17 PC
% updated copyright
%
%6/24/17 PC
% added output varible regarding whether toolbox is current
%
%6/25/17 PC
% updated warning about using a non-stable release
% had the comparisons switch, so it thought old versions were newer
%
%8/13/17 PC
% fixed typo in error message
%
%8/24/17 PC
% never defined output variable for when there is no internet connection.
%  It will now be empty
%
%6/22/18 PC
% it appears that the html code changed at some point on github, so this
%  script has not been working for a while now... I should thing of a better
%  way to do this
%
%8/5/19 PC
% script no longer relies on pulling info from github. this should work
%  better long term

try
    %pull webpage from github
    urlstr = 'http://peterclayson.com/wp-content/uploads/era_currentversion.txt';
    webraw = webread(urlstr,'text','html');

    %trim the whitespace
    ver = strtrim(webraw);
    
    %format version of the currently installed toolbox
    used_parts = sscanf(eraver,'%d.%d.%d')';

    %format version of the toolbox found online
    found_parts = sscanf(ver,'%d.%d.%d')';

    %compare cmdstan versions
    for ii = 1:3
        if used_parts(ii) < found_parts(ii)
            erainstalledver = 0;
            break;
        elseif used_parts(ii) > found_parts(ii)
            erainstalledver = 2;
            break;
        elseif used_parts(ii) == found_parts(ii)
            erainstalledver = 1;
        end
    end
    
    %let user know what was found
    switch erainstalledver
        case 1
            str = 'You are running the most up-to-date version of the toolbox';
            fprintf('\n%s\n',str);
        case 0
            str = 'There is a new version of the toolbox available on ';
            str = [str... 
                '<a href="matlab:web(''https://github.com/peclayson/ERA_Toolbox/releases/latest'',''-browser'')">Github</a>'];
            fprintf('\n%s\n',str);
        case 2
            warning('You are running a non-stable release of the toolbox');
            fprintf(strcat(...
                'You likely cloned github, rather than installed the latest stable release\n',...
                'As a result, complete functionality cannot be guaranteed\n',...
                'The lastest stable release can be downloaded at'));
            str = '<a href="matlab:web(''https://github.com/peclayson/ERA_Toolbox/releases/latest'',''-browser'')">Github</a>';
            fprintf('\n\n%s\n',str);
    end
catch
    fprintf('\nUnable to connect to internet to check for new releases\n');
    erainstalledver = [];
end

end
