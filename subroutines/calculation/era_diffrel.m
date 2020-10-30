function diffscore = era_diffrel(varargin)
%Calculate difference score reliability of single session data
%
%[idtable] = era_diffrel('bp',id_varcov,'bt',trl_varcov,...
%  'er_var',b_sigma,'obs',[20 40],'CI',.95,'est','dep')
%
%Last Modified 9/8/20
%
%Inputs
% bp - between-person (co)variance estimates from CmdStan (id_varcov)
% bt - between-trial (co)variance estimates from CmdStan (id_varcov)
% er_var - error variance estiamtes from CmdStan (b_sigma)
% obs - observations for measure 1 and 2 [trls trls]
% est - type of estimate to be outputted
%  'dep' - dependability, 'gen' - generalizability
% CI - size of the credible interval in decimal format: .95 = 95%
%
%Outputs
% diffscore - array, which contains
%  ll - lower limit of credible interval for dependability/generalizability
%  pt - point estimate of credible interval for dependability/generalizability
%  ul - upper limit of credible interval for dependability/generalizability
%  bp_cov_ll - lower limit of credible interval for between-trial covariance
%  bp_cov_pt - point estimate of credible interval for between-trial covariance
%  bp_cov_ul - upper limit of credible interval for between-trial covariance
%  bt_cov_ll - lower limit of credible interval for between-trial covariance
%  bt_cov_pt - point estimate of credible interval for between-trial covariance
%  bt_cov_ul - upper limit of credible interval for between-trial covariance
%  wp_cov - within-person (residual) covariance
%  icc_ll - lower limit of credible interval for ICC
%  icc_pt - point estimate of credible interval for ICC
%  icc_ul - upper limit of credible interval for ICC
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
% by Peter Clayson (9/8/20)
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
            'See help era_diffrel for more information about inputs'));
    end
    
    %check if bp was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bp',varargin),1);
    if ~isempty(ind)
        bp = varargin{ind+1};
    else
        error('varargin:bp',... %Error code and associated error
            strcat('WARNING: Between-person (co)variances not specified \n\n',...
            'Please input bp. See help era_diffrel for more information \n'));
    end
    
    %check if bt was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('bt',varargin),1);
    if ~isempty(ind)
        bt = varargin{ind+1};
    else
        error('varargin:bt',... %Error code and associated error
            strcat('WARNING: Between-trial (co)variances not specified \n\n',...
            'Please input bt. See help era_diffrel for more information \n'));
    end
    
    %check if er_var was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('er_var',varargin),1);
    if ~isempty(ind)
        er_var = varargin{ind+1};
    else
        error('varargin:er_var',... %Error code and associated error
            strcat('WARNING: Error variance not specified \n\n',...
            'Please input er_var. See help era_diffrel for more information \n'));
    end
    
    %check if est was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('obs',varargin),1);
    if ~isempty(ind)
        obs = varargin{ind+1};
    else
        error('varargin:obs',... %Error code and associated error
            strcat('WARNING: please input number of trials for each measure \n\n',...
            'Please input obs. See help era_diffrel for more information \n'));
    end
    
    %check if est was specified.
    %If it is not found, set display error.
    ind = find(strcmpi('est',varargin),1);
    if ~isempty(ind)
        est = varargin{ind+1};
    else
        error('varargin:est',... %Error code and associated error
            strcat('WARNING: reliability estimate not specified \n\n',...
            'Please input est. See help era_diffrel for more information \n'));
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
                'See help era_diffrel for more information \n'));
        end
    else
        error('varargin:ci',... %Error code and associated error
            strcat('WARNING: Size of credible interval not specified \n\n',...
            'Please input CI. See help era_diffrel for more information \n'));
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

%parse (co)variance components
bp1 = bp(:,1,1);
bp2 = bp(:,2,2);
bp_cov = bp(:,2,1);

bt1 = bt(:,1,1);
bt2 = bt(:,2,2);
bt_cov = bt(:,2,1);

wp1 = exp(er_var(:,1)) .^ 2;
wp2 = exp(er_var(:,2)) .^ 2;

obs1 = obs(1);
obs2 = obs(2);

uni = bp1 + bp2 - (2*bp_cov);

rel_err = (wp1 ./ obs1) + (wp2 ./ obs2);
abs_err = rel_err + (bt1 ./ obs1) + (bt2 ./ obs2) -...
    ((2 .* bt_cov) ./ harmmean(obs));

denom_rel_err = wp1 + wp2;
denom_abs_err = denom_rel_err + bt1 + bt2 - (2 .* bt_cov);

switch est
    case 'dep'
        
        temp = quantile(uni ./ (uni + abs_err),[ciedge, 1-ciedge]);
        ll = temp(1);
        ul = temp(2);
        pt = mean(uni ./ (uni + abs_err));
        
        temp = quantile(uni ./ (uni + denom_abs_err),[ciedge 1-ciedge]);
        icc_ll = temp(1);
        icc_ul = temp(2);
        icc_pt = mean(uni ./ (uni + denom_abs_err));
        
        est_type = 'dependability';
        
    case 'gen'
        
        temp = quantile(uni ./ (uni + rel_err),[ciedge 1-ciedge]);
        ll = temp(1);
        ul = temp(2);
        pt = mean(uni ./ (uni + rel_err));
        
        temp = quantile(uni ./ (uni + denom_rel_err),[ciedge 1-ciedge]);
        icc_ll = temp(1);
        icc_ul = temp(2);
        icc_pt = mean(uni ./ (uni + denom_rel_err));
        
        est_type = 'generalizability';
        
end


%for now, covariance is not estimated for within-person error because the
%two events are not considered to be jointly observed


%Outputs
% diffscore - array, which contains
%  ll - lower limit of credible interval for dependability/generalizability
%  pt - point estimate of credible interval for dependability/generalizability
%  ul - upper limit of credible interval for dependability/generalizability
%  bp_cov_ll - lower limit of credible interval for between-trial covariance
%  bp_cov_pt - point estimate of credible interval for between-trial covariance
%  bp_cov_ul - upper limit of credible interval for between-trial covariance
%  bt_cov_ll - lower limit of credible interval for between-trial covariance
%  bt_cov_pt - point estimate of credible interval for between-trial covariance
%  bt_cov_ul - upper limit of credible interval for between-trial covariance
%  wp_cov - within-person (residual) covariance
%  icc_ll - lower limit of credible interval for ICC
%  icc_pt - point estimate of credible interval for ICC
%  icc_ul - upper limit of credible interval for ICC
% 

%create output structure array
diffscore = [];

diffscore.est_type = est_type;

diffscore.pt = pt;
diffscore.ll = ll;
diffscore.ul = ul;

diffscore.bp_cov_pt = mean(bp_cov);
diffscore.bp_cov_ll = quantile(bp_cov,ciedge);
diffscore.bp_cov_ul = quantile(bp_cov,1-ciedge);

%output between-trial variances as well
diffscore.bt1_var_pt = mean(bt1);
diffscore.bt1_var_ll = quantile(bt1,ciedge);
diffscore.bt1_var_ul = quantile(bt1,1-ciedge);

diffscore.bt2_var_pt = mean(bt2);
diffscore.bt2_var_ll = quantile(bt2,ciedge);
diffscore.bt2_var_ul = quantile(bt2,1-ciedge);

diffscore.bt_cov_pt = mean(bt_cov);
diffscore.bt_cov_ll = quantile(bt_cov,ciedge);
diffscore.bt_cov_ul = quantile(bt_cov,1-ciedge);

%two observations are not considered to be jointly observed
diffscore.wp_cov = 0;

diffscore.icc_pt = icc_pt;
diffscore.icc_ll = icc_ll;
diffscore.icc_ul = icc_ul;

end