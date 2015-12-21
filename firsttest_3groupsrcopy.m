
filename = '/Users/Peter/Documents/RelOutpus/ern_cnt.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


rawNumericColumns = raw(:, [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]);
rawCellColumns = raw(:, 5);

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

erncnt = table;
erncnt.VarName1 = cell2mat(rawNumericColumns(:, 1));
erncnt.id = cell2mat(rawNumericColumns(:, 2));
erncnt.trial = cell2mat(rawNumericColumns(:, 3));
erncnt.urtrial = cell2mat(rawNumericColumns(:, 4));
erncnt.error_label = rawCellColumns(:, 1);
erncnt.error_code = cell2mat(rawNumericColumns(:, 5));
erncnt.ern_fcz = cell2mat(rawNumericColumns(:, 6));
erncnt.ern_roi = cell2mat(rawNumericColumns(:, 7));
erncnt.mdd = cell2mat(rawNumericColumns(:, 8));
erncnt.gad = cell2mat(rawNumericColumns(:, 9));
erncnt.ocd = cell2mat(rawNumericColumns(:, 10));
erncnt.pd = cell2mat(rawNumericColumns(:, 11));
erncnt.ptsd = cell2mat(rawNumericColumns(:, 12));
erncnt.sp = cell2mat(rawNumericColumns(:, 13));
erncnt.dep = cell2mat(rawNumericColumns(:, 14));
erncnt.anx = cell2mat(rawNumericColumns(:, 15));
erncnt.cnt = cell2mat(rawNumericColumns(:, 16));
erncnt.drop = cell2mat(rawNumericColumns(:, 17));
erncnt.group = cell2mat(rawNumericColumns(:, 18));
erncnt.half = cell2mat(rawNumericColumns(:, 19));
erncnt.old_id = cell2mat(rawNumericColumns(:, 20));
erncnt.id2 = cell2mat(rawNumericColumns(:, 21));

clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;


filename = '/Users/Peter/Documents/RelOutpus/ern_anx.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


rawNumericColumns = raw(:, [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]);
rawCellColumns = raw(:, 5);

ernanx = table;
ernanx.VarName1 = cell2mat(rawNumericColumns(:, 1));
ernanx.id = cell2mat(rawNumericColumns(:, 2));
ernanx.trial = cell2mat(rawNumericColumns(:, 3));
ernanx.urtrial = cell2mat(rawNumericColumns(:, 4));
ernanx.error_label = rawCellColumns(:, 1);
ernanx.error_code = cell2mat(rawNumericColumns(:, 5));
ernanx.ern_fcz = cell2mat(rawNumericColumns(:, 6));
ernanx.ern_roi = cell2mat(rawNumericColumns(:, 7));
ernanx.mdd = cell2mat(rawNumericColumns(:, 8));
ernanx.gad = cell2mat(rawNumericColumns(:, 9));
ernanx.ocd = cell2mat(rawNumericColumns(:, 10));
ernanx.pd = cell2mat(rawNumericColumns(:, 11));
ernanx.ptsd = cell2mat(rawNumericColumns(:, 12));
ernanx.sp = cell2mat(rawNumericColumns(:, 13));
ernanx.dep = cell2mat(rawNumericColumns(:, 14));
ernanx.anx = cell2mat(rawNumericColumns(:, 15));
ernanx.cnt = cell2mat(rawNumericColumns(:, 16));
ernanx.drop = cell2mat(rawNumericColumns(:, 17));
ernanx.group = cell2mat(rawNumericColumns(:, 18));
ernanx.half = cell2mat(rawNumericColumns(:, 19));
ernanx.old_id = cell2mat(rawNumericColumns(:, 20));
ernanx.id2 = cell2mat(rawNumericColumns(:, 21));

clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;


filename = '/Users/Peter/Documents/RelOutpus/ern_dep.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

rawNumericColumns = raw(:, [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]);
rawCellColumns = raw(:, 5);


erndep = table;
erndep.VarName1 = cell2mat(rawNumericColumns(:, 1));
erndep.id = cell2mat(rawNumericColumns(:, 2));
erndep.trial = cell2mat(rawNumericColumns(:, 3));
erndep.urtrial = cell2mat(rawNumericColumns(:, 4));
erndep.error_label = rawCellColumns(:, 1);
erndep.error_code = cell2mat(rawNumericColumns(:, 5));
erndep.ern_fcz = cell2mat(rawNumericColumns(:, 6));
erndep.ern_roi = cell2mat(rawNumericColumns(:, 7));
erndep.mdd = cell2mat(rawNumericColumns(:, 8));
erndep.gad = cell2mat(rawNumericColumns(:, 9));
erndep.ocd = cell2mat(rawNumericColumns(:, 10));
erndep.pd = cell2mat(rawNumericColumns(:, 11));
erndep.ptsd = cell2mat(rawNumericColumns(:, 12));
erndep.sp = cell2mat(rawNumericColumns(:, 13));
erndep.dep = cell2mat(rawNumericColumns(:, 14));
erndep.anx = cell2mat(rawNumericColumns(:, 15));
erndep.cnt = cell2mat(rawNumericColumns(:, 16));
erndep.drop = cell2mat(rawNumericColumns(:, 17));
erndep.group = cell2mat(rawNumericColumns(:, 18));
erndep.half = cell2mat(rawNumericColumns(:, 19));
erndep.old_id = cell2mat(rawNumericColumns(:, 20));
erndep.id2 = cell2mat(rawNumericColumns(:, 21));

clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;

% 
% %prepare data
% erncontrols = sortrows(erncontrols,'id');
% id2 = zeros(0,length(erncontrols.id));
% 
% for i = 1:length(erncontrols.id)
%  
%     if i == 1
%         id2(1) = 1;
%         count = 1;
%     elseif i > 1 && strcmp(char(erncontrols.id(i)), char(erncontrols.id(i-1)))
%         id2(i) = count;
%     elseif i > 1 && ~strcmp(erncontrols.id(i), erncontrols.id(i-1))
%         count = count+1;
%         id2(i) = count;
%     end
%     
% end
% 
% erncontrols.id2 = id2(:);


niter = 100;
nchains = 3;

ern_stan = {
  'data {' 
  '  int<lower=0> NDEP; //number of obs in dep'
  '  int<lower=0> NANX; //number of obs in anx'
  '  int<lower=0> NCNT; //number of obs in cnt'
  '  int<lower=0> JDEP; //number of subj in dep'
  '  int<lower=0> JANX; //number of subj in anx'
  '  int<lower=0> JCNT; //number of subj in cnt'
  '  int<lower=0, upper=JDEP> id_dep[NDEP]; //subject id dep;'
  '  int<lower=0, upper=JANX> id_anx[NANX]; //subject id anx;'
  '  int<lower=0, upper=JCNT> id_cnt[NCNT]; //subject id cnt;'
  '  real roi_dep[NDEP];'
  '  real roi_anx[NANX];'
  '  real roi_cnt[NCNT];'
  '}'
  'parameters {'
  '  real mu_dep;'
  '  real mu_anx;'
  '  real mu_cnt;'
  '  real<lower=0> sig_u_dep;'
  '  real<lower=0> sig_u_anx;'
  '  real<lower=0> sig_u_cnt;'
  '  real<lower=0> sig_e_dep;'
  '  real<lower=0> sig_e_anx;'
  '  real<lower=0> sig_e_cnt;'
  '  vector[JDEP] u_raw_dep;'
  '  vector[JANX] u_raw_anx;'
  '  vector[JCNT] u_raw_cnt;'
  '}'
  'transformed parameters {'
  '  vector[JDEP] u_dep;'
  '  vector[JANX] u_anx;'
  '  vector[JCNT] u_cnt;'
  '  u_dep <- mu_dep + sig_u_dep*u_raw_dep;'
  '  u_anx <- mu_anx + sig_u_anx*u_raw_anx;'
  '  u_cnt <- mu_cnt + sig_u_cnt*u_raw_cnt;'
  '}'
  'model {'
  '  u_raw_dep ~ normal(0,1);'
  '  for (i in 1:NDEP) {'
  '    roi_dep[i] ~ normal(u_dep[id_dep[i]], sig_e_dep);'
  '  }'
  '  u_raw_anx ~ normal(0,1);'
  '  for (i in 1:NANX) {'
  '    roi_anx[i] ~ normal(u_anx[id_anx[i]], sig_e_anx);'
  '  }'
  '  u_raw_cnt ~ normal(0,1);'
  '  for (i in 1:NCNT) {'
  '    roi_cnt[i] ~ normal(u_cnt[id_cnt[i]], sig_e_cnt);'
  '  }'
  '  '
  '  mu_dep ~ normal(0,100);'
  '  mu_anx ~ normal(0,100);'
  '  mu_cnt ~ normal(0,100);'
  '  sig_u_dep ~ cauchy(0,40);'
  '  sig_u_anx ~ cauchy(0,40);'
  '  sig_u_cnt ~ cauchy(0,40);'
  '  sig_e_dep ~ cauchy(0,40);'
  '  sig_e_anx ~ cauchy(0,40);'
  '  sig_e_cnt ~ cauchy(0,40);'
  '}'
};


data = struct('NDEP',length(erndep.id), ... %number of observations
    'NANX',length(ernanx.id), ... %number of observations
    'NCNT',length(erncnt.id), ... %number of observations
    'JDEP', length(unique(erndep.id)),... %number of participants
    'JANX', length(unique(ernanx.id)),... %number of participants
    'JCNT', length(unique(erncnt.id)),... %number of participants
    'id_dep', erndep.id2,... %id variable
    'id_anx', ernanx.id2,... %id variable
    'id_cnt', erncnt.id2,... %id variable
    'roi_dep', erndep.ern_roi,...
    'roi_anx', ernanx.ern_roi,...
    'roi_cnt', erncnt.ern_roi); %measurement variable
 
% ern_stan_parms <- c('mu_cnt',
%                     'sig_u_cnt',
%                     'sig_e_cnt')

fit = stan('model_code', ern_stan, 'model_name','test1', 'data', data, 'iter', niter, ...
    'chains', nchains, 'refresh', niter/10,'verbose',true,...
    'file_overwrite', true);

print(fit);



mu = fit.extract('pars','mu_dep').mu_dep;
sig_u = fit.extract('pars','sig_u_dep').sig_u_dep;
sig_e = fit.extract('pars','sig_e_dep').sig_e_dep;
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
hline = refline([0,.7]);
hline.Color = 'r';


