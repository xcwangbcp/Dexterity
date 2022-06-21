clear;close all
monkey_name = {'2','76','132','43','137','187'};
monkey_num  = length(monkey_name);
% RT_mean     = zeros(monkey_num,1);
% RT_std      = RT_mean;
% errate_mean = RT_mean;
% errate_std  = RT_mean;
% control={'95','2','76','35'};nose={'133','132','13'};stritum={'187','43','44'};
% gastrointestional={'137','195','159','25','60','70'};
ctrl_RT =[];  
nose_RT =[];
stri_RT =[];
gas_RT  =[];

crtl_RT_std = [];
nose_RT_std = [];
stri_RT_std = [];
gas_RT_std  = [];

% errate_mean
for i=1:monkey_num
    name = monkey_name(i);
    [RT_mean,RT_std,errate_mean,errate_std] = singleMonkey(name);
    switch char(name)
        case {'2','76','35','95'}
           ctrl_RT     = [ctrl_RT,RT_mean]; 
           crtl_RT_std = [crtl_RT_std,RT_std];
        case {'133','132','13'}
           nose_RT     = [nose_RT,RT_mean];
           nose_RT_std = [nose_RT_std,RT_std];
        case {'137','195','159','25','60','70'}
           gas_RT     = [gas_RT,RT_mean];
           gas_RT_std = [gas_RT_std,RT_std];
        case {'187','43','44'}
           stri_RT     = [stri_RT,RT_mean];
           stri_RT_std = [stri_RT_std,RT_std];
    end
end
RT     = [mean(ctrl_RT),mean(nose_RT),mean(gas_RT),mean(stri_RT)];
RT_std = [mean(crtl_RT_std),mean(nose_RT_std),mean(stri_RT_std),mean(gas_RT_std)];
figure
bar(1:4,RT)
hold on
errorbar(1:4, RT,RT_std)
groupname = {'Control','Nose','Stritum','Gastroint'};
figure
subplot(211)
errorbar(RT,RT_std,'r')
% axis([0 6 200 1000])
title('reaction time in ms')
xticklabels(groupname)
subplot(212)
errorbar(errate_mean,errate_std,'b')
xticks([1 2 3 4 5])
xticklabels(monkey_name)
axis([0 6 0 1.2])
title('correct rate %')


function [RT_mean,RT_std,erro_rat,errate_std] = singleMonkey(monkey_name) 
list        = ls([char(monkey_name) '-*.mat']);
length_file = size(list);
length_file = length_file(1);
erro_rate   = zeros(length_file,1);
RT          = [];
for i= 1:length_file
    R            = load(list(i,:));
    erro_rate(i) = R.R.erro_rate;
    R.R.RT_all   = R.R.RT_all(R.R.RT_all<800);
    RT           = [RT;R.R.RT_all];
end 
erro_rat  = mean(erro_rate);
errate_std= std(erro_rate); 
RT_mean   = mean(RT);
RT_std    = std(RT);
% figure
% errorbar( RT_mean,RT_std,'xr')
end