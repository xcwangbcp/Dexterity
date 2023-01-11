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
dataType   ='normal';
failIndex  = [];
speedIndex = [];
dropIndex  = [];

for yearNum=2:2 %length(yearId)
    k=1;l=1;m=1;n=1;
    year = yearId(yearNum);
    if yearNum==1
        monthId = {'06','07','08','09','10','11','12'};
    else
        monthId = {'01','02','03','04','05','06','07','08','09'};
    end
    monthNum=length(monthId);
    for i=1:monthNum
        month = monthId(i);
        for j = 1: length(monkeyId)
            fileName = cell2mat([monkeyId(j) '-' month '-' year '.xlsx']);
            name = monkeyId(j);
            if exist(fileName)
                [speed,fail,drop] = singleMonth(fileName,dataType);
                switch char(name)
                    case {'02','76','35','95','11'} % control

                        ctrl.failIndex(k)  = [failIndex,fail];
                        ctrl.speedIndex(k) = [speedIndex,speed];
                        ctrl.dropIndex(k)  = [dropIndex,drop];
                        k=k+1;
                    case {'133','132','13'}   % nose
                        nose.failIndex(l)  = [failIndex,fail];
                        nose.speedIndex(l) = [speedIndex,speed];
                        nose.dropIndex(l)  = [dropIndex,drop]; 
                        l=l+1;
                    case {'137','195','159','25','60','70'}%gas
                        gas.failIndex(m)  = [failIndex,fail];
                        gas.speedIndex(m) = [speedIndex,speed];
                        gas.dropIndex(m)  = [dropIndex,drop];        
                        m=m+1;
                    case {'187','43','44'}      %stritum
                        stri.failIndex(n)  = [failIndex,fail];
                        stri.speedIndex(n) = [speedIndex,speed];
                        stri.dropIndex(n)  = [dropIndex,drop];
                        n=n+1;
                end
            end 
        end
    end
%    data.ctrl(yearNum,:)=ctrl;data.nose(yearNum,:)=nose;data.gas(yearNum,:)=gas;data.stri(yearNum,:)=stri;
   data.ctrl=ctrl;data.nose=nose;data.gas=gas;data.stri=stri;
   if strcmpi(dataType,'raw')
       savefile=cell2mat([year '-raw.mat']);
   else
       savefile=cell2mat([year '-normal.mat']);
   end
%    ctrlSave=round([ctrl.failIndex',ctrl.speedIndex',ctrl.dropIndex'],2);
   writematrix(round([ctrl.failIndex',ctrl.speedIndex',ctrl.dropIndex'],4),[savefile(1:end-4) '.xls'],'Sheet',1);
   writematrix(round([nose.failIndex',nose.speedIndex',nose.dropIndex'],4),[savefile(1:end-4) '.xls'],'Sheet',2);
   writematrix(round([gas.failIndex',gas.speedIndex',gas.dropIndex'],4),   [savefile(1:end-4) '.xls'],'Sheet',3);
   writematrix(round([stri.failIndex',stri.speedIndex',stri.dropIndex'],4),[savefile(1:end-4) '.xls'],'Sheet',4);
   save(savefile,'ctrl','nose','gas','stri')
%    draw(data)
end
function [speedIndex,failIndex,dropIndex] = singleMonth(fileName,dataType)
%     list      = ls('*-05-2021.xlsx');
%     fileName  = ('monkeyID','month','year');
%     [fileL,~] = size(list);
%     for fileNum=1:fileL
        tableMonkey  = xlsread(fileName);
        distance     = tableMonkey(5,1);
        slitRate     = tableMonkey(1,1);
        wandRate     = tableMonkey(2,1);
        dropRate     = tableMonkey(3,1);
        totalRate    = (slitRate+wandRate);
        if dropRate>0.35
           dropRate=0.35;
        end
        
        if  totalRate>1
            totalRate=1;
        end
        if strcmpi(dataType,'raw')
            fetchtT      = tableMonkey(4,1);
            dropIndex    = dropRate;%/distance;
            speedIndex   = fetchtT;%/distance;
            failIndex    = totalRate;%/distance;
        else
            fetchtT      = tableMonkey(4,1)/distance;
            dropIndex    = dropRate/distance;
            speedIndex   = fetchtT/distance;
            failIndex    = totalRate/distance;
        end

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
%     hold on
    subplot(311)
    bar(1:4,speedIndex,0.4)%'FaceColor',[1 1 1])
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
end
function plotf(data)
clear;
data21=load('2021raw.mat'); data22=load('2022raw.mat');
ctrl21 = data21.ctrl; nose21 = data21.nose; gas21 = data21.gas; stri21 = data21.stri;
ctrl22 = data22.ctrl; nose22 = data22.nose; gas22 = data22.gas; stri22 = data22.stri;
    speedIndex  = [mean(ctrl21.speedIndex) mean(ctrl22.speedIndex); mean(nose21.speedIndex) mean(nose22.speedIndex); ...
                  mean(gas21.speedIndex) mean(gas22.speedIndex); mean(stri21.speedIndex) mean(stri22.speedIndex)];
    
    speed_std   = [std(ctrl21.speedIndex) std(ctrl22.speedIndex);std(nose21.speedIndex) std(nose22.speedIndex);...
                   std(gas21.speedIndex) std(gas22.speedIndex);std(stri21.speedIndex) std(stri22.speedIndex)];

    failIndex   = [mean(ctrl21.failIndex) mean(ctrl22.failIndex); mean(nose21.failIndex) mean(nose22.failIndex); ...
                  mean(gas21.failIndex) mean(gas22.failIndex); mean(stri21.failIndex) mean(stri22.failIndex)];

    fail_std    = [std(ctrl21.failIndex) std(ctrl22.failIndex);std(nose21.failIndex) std(nose22.failIndex);...
                  std(gas21.failIndex) std(gas22.failIndex);std(stri21.failIndex) std(stri22.failIndex)];

    dropIndex   = [mean(ctrl21.dropIndex) mean(ctrl22.dropIndex); mean(nose21.dropIndex) mean(nose22.dropIndex); ...
                  mean(gas21.dropIndex) mean(gas22.dropIndex); mean(stri21.dropIndex) mean(stri22.dropIndex)];

    drop_std    = [std(ctrl21.dropIndex) std(ctrl22.dropIndex);std(nose21.dropIndex) std(nose22.dropIndex);...
                  std(gas21.dropIndex) std(gas22.dropIndex);std(stri21.dropIndex) std(stri22.dropIndex)];
     figure
%     hold on
    positon=[0.85 1.15; 1.85 2.15; 2.85 3.15;3.85 4.15];
    subplot(311)
    bar(speedIndex,0.4)%,'FaceColor',[1 1 1])
    hold on
    errorbar(positon,speedIndex,speed_std,'o','linewidth',1.5,'LineStyle','none')
    groupname = {'Control','Nose','Gastroint','Stritum'};
    xticklabels(groupname)
%     axis([0 5  0 100 ])
    title('average fetch time ')
    ylabel('ms')
    
    subplot(312)
    bar(1:4,failIndex,0.4)%,'FaceColor',[1 1 1])
    hold on
    errorbar(positon, failIndex,fail_std,'o','linewidth',1.5,'LineStyle','none')
    xticklabels(groupname)
%     axis([0 5 0 0.1])
    title('average fail rate')
    ylabel('%')
%     xlabel('Group behavior data on apple-Machine in May 26/28/31')

    subplot(313)
    bar(1:4,dropIndex,0.4)%,'FaceColor',[1 1 1])
    hold on
    errorbar(positon, dropIndex,drop_std,'o','linewidth',1.5,'LineStyle','none')
    xticklabels(groupname)
%     axis([0 5 0 0.015])
    title('average drop rate')
    ylabel('%')
%     xlabel('Group behavior data on apple-Machine in May 26/28/31')  

end