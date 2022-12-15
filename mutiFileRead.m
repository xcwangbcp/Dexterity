% read the single monkey's data and plot one parameter 
% list            = ls(['2' '-*.xlsx']);% list all the files with #2 monkeys
% testDate        = list(:,3:end-5);
% A=datetime(testDate,'Format','MM-yyyy');
% [~,sortDate]    = sort(A);
% fileName        = list(sortDate,:);
% l=size(fileName);
% for fileNum = 1:l(1)
%     table   = xlsread(fileName(fileNum,:));
%     
%     plot(fileNum,table(5,:),'marker','*')
%     hold on
% end
% monthNum21  = ls('*-2021');
% monthNum22   =ls('*-2022');
% monthNum  =  monthNum21+monthNum22;
% monthID   =;
% monthNum  = 14;
clear
% data       = cell(1,4);% 4 groups
yearId     = {'2021','2022'};     
monkeyId   = {'02','25','35','43','44','60','70','76','132','133','137','159','187','195'};
failIndex  = [];
speedIndex = [];
dropIndex  = [];
k=1;l=1;m=1;n=1;
for yearNum=2:2 %length(yearId)
    year = yearId(yearNum);
    if yearNum==1
        monthId = {'05','06','07','08','09','10','11','12'};
    else
        monthId = {'01','02','03','04','05','06'};
    end
    monthNum=length(monthId);
    for i=1:monthNum
        month = monthId(i);
        for j = 1: length(monkeyId)
            fileName = cell2mat([monkeyId(j) '-' month '-' year '.xlsx']);
            name = monkeyId(j);
            if exist(fileName)
                [speed,fail,drop] = singleMonth(fileName);
                switch char(name)
                    case {'02','76','35','95'} % control

                        ctrl.failIndex(yearNum,k)  = [failIndex,fail];
                        ctrl.speedIndex(yearNum,k) = [speedIndex,speed];
                        ctrl.dropIndex(yearNum,k)  = [dropIndex,drop];
                        k=k+1;
                    case {'133','132','13'}   % nose
                        nose.failIndex(yearNum,l)  = [failIndex,fail];
                        nose.speedIndex(yearNum,l) = [speedIndex,speed];
                        nose.dropIndex(yearNum,l)  = [dropIndex,drop]; 
                        l=l+1;
                    case {'137','195','159','25','60','70'}%gas
                        gas.failIndex(yearNum,m)  = [failIndex,fail];
                        gas.speedIndex(yearNum,m) = [speedIndex,speed];
                        gas.dropIndex(yearNum,m)  = [dropIndex,drop];        
                        m=m+1;
                    case {'187','43','44'}      %stritum
                        stri.failIndex(yearNum,n)  = [failIndex,fail];
                        stri.speedIndex(yearNum,n) = [speedIndex,speed];
                        stri.dropIndex(yearNum,n)  = [dropIndex,drop];
                        n=n+1;
                end
            end 
        end
    end
   data.ctrl=ctrl;data.nose=nose;data.gas=gas;data.stri=stri;
   draw(data)
end
function [speedIndex,failIndex,dropIndex] = singleMonth(fileName)
%     list      = ls('*-05-2021.xlsx');
%     fileName  = ('monkeyID','month','year');
%     [fileL,~] = size(list);
%     for fileNum=1:fileL
        tableMonkey  = xlsread(fileName);
        distance     = tableMonkey(5,1);
        slitRate     = tableMonkey(1,1);
        wandRate     = tableMonkey(2,1);
        dropRate     = tableMonkey(3,1);
        fetchtT      = tableMonkey(4,1);
        dropIndex    = 100*dropRate/distance;
        speedIndex   = fetchtT/distance;
        failIndex    = 100*(slitRate+wandRate)/distance;
%     end
end
function draw(data)
    ctrl = data.ctrl; nose = data.nose; gas = data.gas; stri = data.stri;
    speedIndex  = [mean(ctrl.speedIndex), mean(nose.speedIndex) ...
                  mean(gas.speedIndex), mean(stri.speedIndex)];
    
    speed_std   = [std(ctrl.speedIndex),std(nose.speedIndex),...
                   std(gas.speedIndex),std(stri.speedIndex)];

    failIndex   = [mean(ctrl.failIndex), mean(nose.failIndex) ...
                  mean(gas.failIndex), mean(stri.failIndex)];

    fail_std    = [std(ctrl.failIndex),std(nose.failIndex),...
                  std(gas.failIndex),std(stri.failIndex)];

    dropIndex   = [mean(ctrl.dropIndex), mean(nose.dropIndex) ...
                  mean(gas.dropIndex), mean(stri.dropIndex)];

    drop_std    = [std(ctrl.dropIndex),std(nose.dropIndex),...
                  std(gas.dropIndex),std(stri.dropIndex)];
    figure
    subplot(311)
    bar(1:4,speedIndex,0.4,'FaceColor',[1 1 1])
    hold on
    errorbar(1:4, speedIndex,speed_std,'o','linewidth',1.5,'LineStyle','none')
    groupname = {'Control','Nose','Gastroint','Stritum'};
    xticklabels(groupname)
    axis([0 5  0 100 ])
    title('speedIndex')
    ylabel('s/m')
    
    subplot(312)
    bar(1:4,failIndex,0.4,'FaceColor',[1 1 1])
    hold on
    errorbar(1:4, failIndex,fail_std,'o','linewidth',1.5,'LineStyle','none')
    xticklabels(groupname)
    axis([0 5 0 10])
    title('failIndex')
    ylabel('%/m')
%     xlabel('Group behavior data on apple-Machine in May 26/28/31')

    subplot(313)
    bar(1:4,dropIndex,0.4,'FaceColor',[1 1 1])
    hold on
    errorbar(1:4, dropIndex,drop_std,'o','linewidth',1.5,'LineStyle','none')
    xticklabels(groupname)
    axis([0 5 0 1])
    title('dropIndex')
    ylabel('%/m')
%     xlabel('Group behavior data on apple-Machine in May 26/28/31')


    savefile='2021-5.';
    save(savefile,'ctrl','nose','gas','stri')
end