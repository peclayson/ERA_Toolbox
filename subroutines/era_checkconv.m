function RELout = era_checkconv(REL)
%Check for convergence by examining the potential scale reduction factor
%(r_hat) and the effective sample size (n_eff)
%
%[REL, rerun] = era_checkconv(REL)
%
%Lasted Updated 1/19/17
%
%Required Input:
% REL - structure array created by era_computerel
%
%Outputs:
% RELout - structure array with the following added fields
%  out.conv.converged - 1 data converged, 0 data did not converge

% Copyright (C) 2016-2018 Peter E. Clayson
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
%

%make sure an input was provided
if length(nargin) ~= 1
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Input not specified \n\n',... 
    'See help era_checkconv for more information on  inputs'));
end

%pull the dimensions to figure out if there are nested cell arrays
%if there are there will be 1 row and multiple columns (a = 1)
[a,nloops] = size(REL.out.conv.data);

%if there is more than 1 row then that means there is only one set of
%convergence statistics to examine
if a == 1

    for i=1:nloops
        %pull data from REL
        %first check r_hats
        data = REL.out.conv.data{i};
        rhats = data(:,[1 3]);
        
        %see if any of the rhats didn't equal 1 (i.e., did not converge)
        indbad = find([rhats{2:end,2}] >= 1.1);

        %specify whether convergenece between chains was reached
        if ~isempty(indbad)
            result = 0;
            break;
        else
            result = 1;
        end

        if result == 1
            %check whether neff is (5*2*nchains)
            neff = data(:,[1 2]);

            %see if any of the rhats didn't equal 1 (i.e., did not converge)
            indbad = find([neff{2:end,2}] < (5*2*REL.nchains));

            %specify whethere convergenece between chains was reached
            if ~isempty(indbad)
                result = 0;
                break;
            else
                result = 1;
            end
        end
    end
    
elseif a ~= 1
    
    %pull data from REL
    %first check r_hats
    rhats = REL.out.conv.data(:,[1 3]);

    %see if any of the rhats didn't equal 1 (i.e., did not converge)
    indbad = find([rhats{2:end,2}] >= 1.1);

    %specify whether convergenece between chains was reached
    if ~isempty(indbad)
        result = 0;
    else
        result = 1;
    end

    if result == 1
        %check whether neff is (5*2*nchains)
        neff = REL.out.conv.data(:,[1 2]);

        %see if any of the rhats didn't equal 1 (i.e., did not converge)
        indbad = find([neff{2:end,2}] < (5*2*REL.nchains));

        %specify whethere convergenece between chains was reached
        if ~isempty(indbad)
            result = 0;
        else
            result = 1;
        end
    end
end
    
    
REL.out.conv.converged = result;
RELout = REL;

end

