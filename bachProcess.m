clear;close all
monkey_name = {'2','76','132','43','137','187'};
monkey_num  = length(monkey_name);
RT_mean     = zeros(monkey_num,1);
RT_std      = RT_mean;
errate_mean = RT_mean;
errate_std  = RT_mean;
for i=1:monkey_num
    name = monkey_name(i);
    [RT_mean(i,1),RT_std(i,1),errate_mean(i),errate_std(i)] = singleMonkey(name);
    switch char(name)
        case 2
        case
    
    end
end
figure
subplot(211)
errorbar(RT_mean,RT_std,'r')
axis([0 6 200 1000])
title('reaction time in ms')
xticklabels(monkey_name)
subplot(212)
errorbar(errate_mean,errate_std,'b')
xticks([1 2 3 4 5])
xticklabels(monkey_name)
axis([0 6 0 1.2])
title('correct rate %')
% control={'95','2','76','35'};nose={'133','132','13'};stritum={'187','43','44'};
% gastrointestional={'137','195','159','25','60','70'};
group_ ={}
figure
bar(1:5,RT_mean)
hold on
errorbar(1:5, RT_mean,RT_std)
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