function ra_relfigures(varargin)





if ~isempty(varargin)
    
    %the optional inputs check assumes that there was an even number of 
    %optional inputs entered. If not, an error will displayed and the
    %script will terminate.
    if mod(length(varargin),2)  
        error('varargin:incomplete',... %Error code and associated error
        strcat('WARNING: Inputs are incomplete \n\n',... 
        'Make sure each variable input is paired with a value \n',...
        'See help ra_relfigures for more information on optional inputs'));
    end
    
    %check if a location for the file to be loaded was specified. 
    %If it is not found, set display error.
    ind = find(strcmp('data',varargin),1);
    if ~isempty(ind)
        REL = varargin{ind+1}; 
    else 
        error('varargin:nofile',... %Error code and associated error
        strcat('WARNING: File location not specified \n\n',... 
        'Please input the full path specifying the file to be loaded \n'));
    end
   
elseif ~isempty(varargin)
    
    error('varargin:incomplete',... %Error code and associated error
    strcat('WARNING: Optional inputs are incomplete \n\n',... 
    'Make sure each variable input is paired with a value \n',...
    'See help ra_relfigures for more information on optional inputs'));
    
end %if ~isempty(varargin)

data = struct;
data.label = {};
data.mu = [];
data.mu.raw = [];
data.mu.ll = [];
data.mu.ul = [];
data.sig_u = [];
data.sig_u.raw = [];
data.sig_u.ll = [];
data.sig_u.ul = [];
data.sig_e = [];
data.sig_e.raw = [];
data.sig_e.ll =[];
data.sig_e.ul = [];

for i=1:length(REL.out.labels)
    
    data(i).label = REL.out.labels(i);
    data(i).mu.raw = REL.out.mu(:,i);
    data(i).sig_u.raw = REL.out.sig_u(:,i);
    data(i).sig_e.raw = REL.out.sig_e(:,i);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitdata = table(mu, sig_u, sig_e, lp);

%wrote my own icc function
mean_icc = mean(icc(sig_u(:), sig_e(:)));
ll_icc = quantile(icc(sig_u(:), sig_e(:)),.025);
ul_icc = quantile(icc(sig_u(:), sig_e(:)),.975);

mean_sig_u = mean(sig_u);
mean_sig_e = mean(sig_e);
mean_mu = mean(mu);

ll_sig_u = quantile(sig_u,.025);
ul_sig_u = quantile(sig_u,.975);
ll_sig_e = quantile(sig_e,.025);
ul_sig_e = quantile(sig_e,.975);
ll_mu = quantile(mu,.025);
ul_mu = quantile(mu,.975);

mrel = zeros(0,50);
llrel = zeros(0,50);
ulrel = zeros(0,50);
x = zeros(0,50);

for i = 1:50
   
   x(i) = i;
   mrel(i) = mean(reliab(sig_u,sig_e,i)); 
   llrel(i) = quantile(reliab(sig_u,sig_e,i),.025);
   ulrel(i) = quantile(reliab(sig_u,sig_e,i),.975);
   
end

plot(x,mrel);



end

function depout = reliab(var_u, var_e, obs)

depout = var_u.^2 ./ (var_u.^2 + (var_e.^2./obs));

end

function iccout = icc(var_u,var_e)

iccout = var_u.^2 ./ (var_u.^2 + var_e.^2);

end
