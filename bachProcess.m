clear;close all
monkey_name = {'2','76','132','43','137','187'};
monkey_num  = length(monkey_name);
% control={'95','2','76','35'};nose={'133','132','13'};stritum={'187','43','44'};
% gastrointestional={'137','195','159','25','60','70'};
ctrl.RT =[];     ctrl.errate = [];
nose.RT =[];     nose.errate = [];
stri.RT =[];     stri.errate = [];
gas.RT  =[];     gas.errate  = [];

ctrl.RT_std = []; ctrl.er_std = [];
nose.RT_std = []; nose.er_std = []; 
stri.RT_std = []; stri.er_std = [];
gas.RT_std  = []; gas.er_std  = [];


% errate_mean
for i=1:monkey_num
    name = monkey_name(i);
    [RT_mean,RT_std,errate,errate_std] = singleMonkey(name);
    switch char(name)
        case {'2','76','35','95'} % control
           ctrl.RT     = [ctrl.RT,RT_mean]; 
           ctrl.RT_std = [ctrl.RT_std,RT_std];
           ctrl.errate = [ctrl.errate,errate];
           ctrl.er_std = [ctrl.er_std,errate_std];  
        case {'133','132','13'}   % nose
           nose.RT     = [nose.RT,RT_mean];
           nose.RT_std = [nose.RT_std,RT_std];
           nose.errate = [nose.errate,errate];
           nose.er_std = [nose.er_std,errate_std];
        case {'137','195','159','25','60','70'}%gas
           gas.RT      = [gas.RT,RT_mean];
           gas.RT_std  = [gas.RT_std,RT_std];
           gas.errate  = [gas.errate,errate];
           gas.er_std  = [gas.er_std,errate_std];
        case {'187','43','44'}      %stritum
           stri.RT     = [stri.RT,RT_mean];
           stri.RT_std = [stri.RT_std,RT_std];
           stri.errate = [stri.errate,errate];
           stri.er_std = [stri.er_std,errate_std];
    end
end
RT       = [mean(ctrl.RT),    mean(nose.RT),    mean(gas.RT),    mean(stri.RT)];
RT_std   = [mean(ctrl.RT_std),mean(nose.RT_std),mean(gas.RT_std),mean(stri.RT_std),];
rate     = [mean(ctrl.errate),mean(nose.errate),mean(gas.errate),mean(stri.errate)];
rate_std = [mean(ctrl.er_std),mean(nose.er_std),mean(gas.er_std),mean(stri.er_std)];

figure
subplot(211)
bar(1:4,RT,0.4,'FaceColor',[1 1 1])
hold on
errorbar(1:4, RT,RT_std,'o','linewidth',1.5,'LineStyle','none')
groupname = {'Control','Nose','Gastroint','Stritum'};
xticklabels(groupname)
axis([0 5 0 1000])
title('Time to fetch an apple')
ylabel('time in ms')

subplot(212)
bar(1:4,rate,0.4,'FaceColor',[1 1 1])
hold on
errorbar(1:4, rate,rate_std,'o','linewidth',1.5,'LineStyle','none')
xticklabels(groupname)
axis([0 5 0 1.1])
title('Error rate')
ylabel('%')
xlabel('Group behavior data on apple-Machine in May 26/28/31')
savefile='2021-5';
save(savefile,'ctrl','nose','gas','stri')




function [RT_mean,RT_std,erro_rat,errate_std] = singleMonkey(monkey_name) 
list        = ls([char(monkey_name) '-*.mat']);
length_file = size(list);
length_file = length_file(1);
erro_rate   = zeros(length_file,1);
RT          = [];
for i= 1:length_file
    R            = load(list(i,:));
    erro_rate(i) = R.R.erro_rate;

    if erro_rate(i)>1
        erro_rate(i)=1;
    end
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