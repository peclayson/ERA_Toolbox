function idtable = era_ssrel(varargin)
%Calculate subject-level dependability and ICCs using variance components  
% from CmdStan
%
%[idtable] = era_dep('bp',gro_sds,'wp_pop',pop_sdlog,'wp_ss',ind_sdlog,...
%  'idtable',idtable,'CI',.95)
%
%Last Modified 8/26/20
%
%Inputs
% bp - between-person variance components from CmdStan (gro_sds)
% wp_pop - population estimates of within-person variance components from 
%  CmdStan (pop_sdlog)
% wp_ss - subject-level estimates of within-person variance components from 
%  CmdStan (ind_sdlog)
% idtable - table with information about original id, new id2 for CmdStan,
%  and the number of trials retained for each id (each row id should line
%  up exactly with each column for wp_ss)
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Outputs
% tableout - a table that is concatenated with idtable input.
%   id - original id
%   id2 - CmdStan id
%   trls - number of retained trials
%   dep_pt - dependability point estimate
%   dep_ll - lower limit of the credible interval specified by CI for
%     dependability
%   dep_ul - upper limit of the credible interval specified by CI for
%     dependability
%   icc_pt - ICC point estimate
%   icc_ll - lower limit of the credible interval specified by CI for
%     ICC
%   icc_ul - upper limit of the credible interval specified by CI for
%     ICC
%   bp_var - population between-person variance estimate (same for all ids)
%   ss_errvar - subject-level error variance estimates
%   pop_errvar - population-level error variance estimate (same for all 
%     ids)
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
% by Peter Clayson (8/26/20)
% peter.clayson@gmail.com
%


%somersault through inputs
if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help era_dep_ssrel for more information about inputs'));
    end
    
    %check if bp was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1}; 
    else 
        error('varargin:bp',... %Error code and associated error
        strcat('WARNING: Between-person variance not specified \n\n',... 
        'Please input bp. See help era_dep_ssrel for more information \n'));
    end
    
    %check if wp_pop was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('wp_pop',varargin),1);
    if ~isempty(ind)
        wp_pop = varargin{ind+1}; 
    else 
        error('varargin:wp_pop',... %Error code and associated error
        strcat('WARNING: Within-person popluation variance not specified \n\n',... 
        'Please input wp_pop. See help era_dep_ssrel for more information \n'));
    end
    
    %check if wp_ss was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('wp_ss',varargin),1);
    if ~isempty(ind)
        wp_ss = varargin{ind+1}; 
    else 
        error('varargin:wp_ss',... %Error code and associated error
        strcat('WARNING: Subject-level within-person variance not specified \n\n',... 
        'Please input wp_ss. See help era_dep_ssrel for more information \n'));
    end
    
    %check if idtable was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('idtable',varargin),1);
    if ~isempty(ind)
        idtable = varargin{ind+1}; 
    else 
        error('varargin:idtable',... %Error code and associated error
        strcat('WARNING: idtable not specified \n\n',... 
        'Please input idtable. See help era_dep_ssrel for more information \n'));
    end
    
    %check if CI was specified. 
    %If it is not found, set display error.
    ind = find(strcmpi('CI',varargin),1);
    if ~isempty(ind)
        ciperc = varargin{ind+1};
        if ciperc > 1 || ciperc < 0
            error('varargin:ci',... %Error code and associated error
                strcat('WARNING: Size of credible interval should ',...
                'be a value between 0 and 1\n',...
                'A value of ',sprintf(' %2.2f',ciperc),...
                ' is invalid\n',...
                'See help era_dep_ssrel for more information \n'));
        end
    else 
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',... 
            'Please input CI. See help era_dep_ssrel for more information \n'));
    end
    
end

%make sure the user understood the the CI input is width not edges, if
%inputted incorrectly, provide a warning, but don't change
if ciperc < .5 
    str = sprintf(' %2.0f%%',100*ciperc);
    warning('ci:width',... %Warning code and associated warning
        strcat('WARNING: Size of credible interval is small \n\n',...
        'User specified a credible interval width of ',str,'%\n',...
        'If this was intended, ignore this warning.\n'));
end

%calculate edges for CI
ciedge = (1-ciperc)/2;

%estimate between-person and subject-level within-person variance
bp = cell2mat(bp) .^ 2;
wp = exp(wp_pop + cell2mat(wp_ss)) .^ 2;

%estimate subject-level dependability
dep_ll = quantile(bp ./ (bp + (wp ./ idtable.trls(:)')),ciedge);
dep_pt = mean(bp ./ (bp + (wp ./ idtable.trls(:)')));
dep_ul = quantile(bp ./ (bp + (wp ./ idtable.trls(:)')),1-ciedge);

%estimate subject-level ICCs
icc_ll = quantile(bp ./ (bp + wp),ciedge);
icc_pt = mean(bp ./ (bp + wp));
icc_ul = quantile(bp ./ (bp + wp),1-ciedge);

%put information into output table
idtable.dep_pt = dep_pt';
idtable.dep_ll = dep_ll';
idtable.dep_ul = dep_ul';

idtable.icc_pt = icc_pt';
idtable.icc_ll = icc_ll';
idtable.icc_ul = icc_ul';

idtable.bp_var = repmat(mean(bp), height(idtable), 1);
idtable.ss_errvar = mean(wp,1)';
idtable.pop_errvar = repmat(mean(exp(wp_pop).^2),...
    height(idtable), 1);


end