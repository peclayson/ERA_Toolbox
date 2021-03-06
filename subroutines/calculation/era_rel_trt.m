function [ll,pt,ul] = era_rel_trt(varargin)
%Calculate reliability using variance components from CmdStan
%
%[ll,pt,ul] = era_rel_trt('gcoeff',1,'reltype',1,...
% 'bp',sig_u,'bo',sig_p,'bt',sig_t,...
% 'txp',sig_txp,'oxp',sig_oxp,'txo',sig_oxt,'err',sig_err,...
% 'obs',[1 50],'CI',.95)
%
%Last Modified 6/21/19
%
%Inputs
% gcoeff - g-theory coefficient to calculate 1 - dependability, 2 -
%  generalizability
% reltype - reliability coefficient to plot/calculate - 1 - coefficent of
%  equivalence (internal consistency), 2 - coefficient of stability (test
%  retest reliability)
% bp - between-person variance components from CmdStan (sig_id)
% bo - between-occasion variance compoentns from CmdStan (sig_occ)
% bt - between-trial variance components from CmdStan (sig_trl)
% txp - trial x person variance components from CmdStan (sig_trlxid)
% oxp - occasion x person variance components from CmdStan (sig_occxid)
% txo - trial x occasion variance components from CmdStan (sig_trlxocc)
% err - residual variance components from CmdStan (sig_err)
% obs - number of observations. Can either be 1 number or range [min max]
%  E.g., The overall reliability estimate will use one number (central
%   tendency). The plot for the number of trials v reliability will use a
%   range of numbers.
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Outputs
% ll - lower limit of the credible interval specified by CI for
%  dependability
% pt - the point estimate for reliabilty
% ul - upper limit of the credible interval specified by CI for
%  dependability

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
% by Peter Clayson (6/21/19)
% peter.clayson@gmail.com
%
%6/21/19 PC
% updated copyright
%
%10/28/19 PC
% updated formulas 

%somersault through inputs

if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)
        error('varargin:incomplete',... %Error code and associated error
            strcat('WARNING: Inputs are incomplete \n\n',...
            'Make sure each variable input is paired with a value \n',...
            'See help era_rel_trt for more information about inputs'));
    end
    
    %check if gcoeff was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('gcoeff',varargin),1);
    if ~isempty(ind)
        gcoeff = varargin{ind+1};
    else
        error('varargin:gcoeff',... %Error code and associated error
            strcat('WARNING: G-theory estimate not specified \n\n',...
            'Please input gcoeff. See help era_rel_trt for more information \n'));
    end
    
    %check if reltype was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('reltype',varargin),1);
    if ~isempty(ind)
        reltype = varargin{ind+1};
    else
        error('varargin:reltype',... %Error code and associated error
            strcat('WARNING: reliability coefficient not specified \n\n',...
            'Please input reltype. See help era_rel_trt for more information \n'));
    end
    
    %check if bp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1};
    else
        error('varargin:bp',... %Error code and associated error
            strcat('WARNING: Between-person variance not specified \n\n',...
            'Please input bp. See help era_rel_trt for more information \n'));
    end
    
    %check if bo was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bo',varargin),1);
    if ~isempty(ind)
        bo = varargin{ind+1};
    else
        error('varargin:bo',... %Error code and associated error
            strcat('WARNING: Between-occasion variance not specified \n\n',...
            'Please input bo. See help era_rel_trt for more information \n'));
    end
    
    %check if bt was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bt',varargin),1);
    if ~isempty(ind)
        bt = varargin{ind+1};
    else
        error('varargin:bt',... %Error code and associated error
            strcat('WARNING: Between-trial variance not specified \n\n',...
            'Please input bt. See help era_rel_trt for more information \n'));
    end
    
    %check if txp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('txp',varargin),1);
    if ~isempty(ind)
        txp = varargin{ind+1};
    else
        error('varargin:txp',... %Error code and associated error
            strcat('WARNING: trialxperson variance not specified \n\n',...
            'Please input txp. See help era_rel_trt for more information \n'));
    end
    
    %check if oxp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('oxp',varargin),1);
    if ~isempty(ind)
        oxp = varargin{ind+1};
    else
        error('varargin:oxp',... %Error code and associated error
            strcat('WARNING: occasionxperson variance not specified \n\n',...
            'Please input oxp. See help era_rel_trt for more information \n'));
    end
    
    %check if txo was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('txo',varargin),1);
    if ~isempty(ind)
        txo = varargin{ind+1};
    else
        error('varargin:txo',... %Error code and associated error
            strcat('WARNING: trialxoccasion variance not specified \n\n',...
            'Please input txo. See help era_rel_trt for more information \n'));
    end
    
    %check if err was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('err',varargin),1);
    if ~isempty(ind)
        err = varargin{ind+1};
    else
        error('varargin:err',... %Error code and associated error
            strcat('WARNING: error variance not specified \n\n',...
            'Please input err. See help era_rel_trt for more information \n'));
    end
    
    %check if obs was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('obs',varargin),1);
    if ~isempty(ind)
        obs = varargin{ind+1};
    else
        error('varargin:obs',... %Error code and associated error
            strcat('WARNING: Number of observations not specified \n\n',...
            'Please input obs. See help era_rel_trt for more information \n'));
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
                'See help era_rel_trt for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help era_rel_trt for more information \n'));
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

%for now, all trt analyses will use nocc = 1;
nocc = 1;

%inputs are all in standard deviation units, but we need variance
%components to calculatre reliability estimates
%convert inputs to variance components to make formulas cleaner
bp = bp.^2;
bo = bo.^2;
bt = bt.^2;
txp = txp.^2;
oxp = oxp.^2;
txo = txo.^2;
err = err.^2;

if length(obs) == 1
    if gcoeff == 1 &&... %dependability
            reltype == 1 %coefficient of equivalence
        
        %universe score
        uni_scor = bp + (oxp ./ nocc);
        
        %absolute error
        err_term = (txp ./ obs) + (err ./ (obs*nocc)) +...
            (bt ./ obs) + (txo ./ (obs*nocc));
        
        %dependability estimate
        rel = uni_scor ./ (uni_scor + err_term);
        
    elseif gcoeff == 1 &&... %dependabilty
            reltype == 2 %coefficient of stability
        
        %universe score
        uni_scor = bp + (txp ./ obs);
        
        %absolute error
        err_term = (oxp ./ nocc) + (err ./ (obs*nocc)) +...
            (bo ./ nocc) + (txo ./ (obs*nocc));
        
        %dependability estimate
        rel = uni_scor ./ (uni_scor + err_term);
        
        
    elseif gcoeff == 2 &&... %generalizability
            reltype == 1 %coefficient of equivalence
        
        %universe score
        uni_scor = bp + (oxp ./ nocc);
        
        %relative error
        err_term = (txp ./ obs) + (err ./ (obs*nocc));
        
        %generalizability estimate
        rel = uni_scor ./ (uni_scor + err_term);
        
    elseif gcoeff == 2 &&... %generalizability
            reltype == 2 %coefficient of stability
        
        %universe score
        uni_scor = bp + (txp ./ obs);
        
        %relative error
        err_term = (oxp ./ nocc) + (err ./ (obs*nocc));
        
        %generalizability estimate
        rel = uni_scor ./ (uni_scor + err_term);
        
    end
    
    %pull the lower limit, point estimate, and upper limit
    ll = quantile(rel,ciedge);
    pt = mean(rel);
    ul = quantile(rel,1-ciedge);
    
elseif length(obs) == 2 %if obs contains two numbers calculate the
    %reliability for the range specified
    
    %create empty matrices for storing values
    ll = zeros(0,obs(2)-obs(1)+1);
    pt = zeros(0,obs(2)-obs(1)+1);
    ul = zeros(0,obs(2)-obs(1)+1);
    
    for ii = obs(1):obs(2)
        if gcoeff == 1 &&... %dependability
                reltype == 1 %coefficient of equivalence
            
            %universe score
            uni_scor = bp + (oxp ./ nocc);
            
            %absolute error
            err_term = (txp ./ ii) + (err ./ (ii*nocc)) +...
                (bt ./ ii) +  (txo ./ (ii*nocc));
            
            %dependability estimate
            rel = uni_scor ./ (uni_scor + err_term);
            
        elseif gcoeff == 1 &&... %dependabilty
                reltype == 2 %coefficient of stability
            
            %universe score
            uni_scor = bp + (txp ./ ii);
            
            %absolute error
            err_term = (oxp ./ nocc) + (err ./ (ii*nocc)) +...
                (bo ./ nocc) + (txo ./ (ii*nocc));
            
            %dependability estimate
            rel = uni_scor ./ (uni_scor + err_term);
            
            
        elseif gcoeff == 2 &&... %generalizability
                reltype == 1 %coefficient of equivalence
            
            %universe score
            uni_scor = bp + (oxp ./ nocc);
            
            %relative error
            err_term = (txp ./ ii) + (err ./ (ii*nocc));
            
            %generalizability estimate
            rel = uni_scor ./ (uni_scor + err_term);
            
        elseif gcoeff == 2 &&... %generalizability
                reltype == 2 %coefficient of stability
            
            %universe score
            uni_scor = bp + (txp ./ ii);
            
            %relative error
            err_term = (oxp ./ nocc) + (err ./ (ii*nocc));
            
            %generalizability estimate
            rel = uni_scor ./ (uni_scor + err_term);
            
        end
        
        %pull the lower limit, point estimate, and upper limit
        ll(ii) = quantile(rel,ciedge);
        pt(ii) = mean(rel);
        ul(ii) = quantile(rel,1-ciedge);
        
    end %for ii = obs(1):obs(2)
   
end %if length(obs) == 1

end