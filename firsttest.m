
filename = '/Users/Peter/Dropbox/Documents/Publications-In Progress/ERP_ReliabilityRev/NewScripts/ern_controls.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%f%s%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

erncontrols = table(dataArray{1:end-1}, 'VariableNames', {'id','event','ern_roi'});

clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%prepare data
erncontrols = sortrows(erncontrols,'id');
id2 = zeros(0,length(erncontrols.id));

for i = 1:length(erncontrols.id)
 
    if i == 1
        id2(1) = 1;
        count = 1;
    elseif i > 1 && strcmp(char(erncontrols.id(i)), char(erncontrols.id(i-1)))
        id2(i) = count;
    elseif i > 1 && ~strcmp(erncontrols.id(i), erncontrols.id(i-1))
        count = count+1;
        id2(i) = count;
    end
    
end

erncontrols.id2 = id2(:);

clear i id2;

niter = 100;
nchains = 3;

ern_stan = {
  'data {' 
  '  int<lower=0> NCNT; //number of obs'
  '  int<lower=0> JCNT; //number of subj'
  '  int<lower=0, upper=JCNT> id_cnt[NCNT]; //subject id cnt;'
  '  real fcz_cnt[NCNT];'
  '}'
  'parameters {'
  '  real mu_cnt;'
  '  real<lower=0> sig_u_cnt;'
  '  real<lower=0> sig_e_cnt;'
  '  vector[JCNT] u_raw_cnt;'
  '}'
  'transformed parameters {'
  '  vector[JCNT] u_cnt;'
  '  u_cnt <- mu_cnt + sig_u_cnt*u_raw_cnt;'
  '}'
  'model {'
  '  u_raw_cnt ~ normal(0,1);'
  '  for (i in 1:NCNT) {'
  '    fcz_cnt[i] ~ normal(u_cnt[id_cnt[i]], sig_e_cnt);'
  '  }'
  '  '
  '  mu_cnt ~ normal(0,100);'
  '  sig_u_cnt ~ cauchy(0,40);'
  '  sig_e_cnt ~ cauchy(0,40);'
  '}'
};


data = struct('NCNT',length(erncontrols.id2), ... %number of observations
  'JCNT', length(unique(erncontrols.id2)),... %number of participants
  'id_cnt', erncontrols.id2,... %id variable
  'fcz_cnt', erncontrols.ern_fcz); %measurement variable
 
fit = stan('model_code', ern_stan, 'model_name','test1', 'data', data, 'iter', niter, ...
    'chains', nchains, 'refresh', niter/10,'verbose',true,...
    'file_overwrite', true);

fit.block();
print(fit);

mu = fit.extract('pars','mu_cnt').mu_cnt;
sig_u = fit.extract('pars','sig_u_cnt').sig_u_cnt;
sig_e = fit.extract('pars','sig_e_cnt').sig_e_cnt;
lp = fit.extract('pars','lp__').lp__;

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



