close all
clear 

[filename_raw_hand,pathname] = uigetfile('*.csv','Pick a csv file to load in ');
[raw,Txt,~]     = xlsread([pathname,filename_raw_hand]);

% filename_raw_hand  =  'D:\Code\Data\Apple\2-trimmed.csv';
filename_raw_apple = 'D:\Code\Data\Apple\2-trimmed-apple.csv';
% Acrylic_Edge        = 'D:\Code\Dexterity\137LS_trimmed_cam1.csv';
% raw_hand  = table2array(readtable([pathname,filename_raw_hand]));
raw_apple = table2array(readtable(filename_raw_apple));
% Edge_raw  = table2array(readtable(Acrylic_Edge));
Edge(:,1)  = mean(raw_apple(:,1));
Edge(:,2)  = mean(raw_apple(:,4));
Edge       = mean(Edge);

% M.base   = findcolum(Txt,'base_x','base_y','base_p');
% M.MCP1   = findcolum(Txt,'MCP1_x','MCP1_y','MCP1_p');
% M.PIP1   = findcolum(Txt,'PIP1_x','PIP1_y','PIP1_p');
% M.tip1   = findcolum(Txt,'tip1_x','tip1_y','tip1_p'); % thumb tip
% M.MCP2   = findcolum(Txt,'MCP2_x','MCP2_y','MCP2_p');
% M.PIP2   = findcolum(Txt,'PIP2_x','PIP2_y','PIP2_p');
% M.DIP2   = findcolum(Txt,'DIP2_x','DIP2_y','DIP2_p');
M.tip2   = findcolum(Txt,'tip2_x','tip2_y','tip2_p'); % indexfinger tip 
index_tip= raw(:,M.tip2(1));
% index_tip= movmean(index_tip,9);
apple_p    = raw_apple(:,end);
apple      = raw_apple(:,7);
apple(apple_p<0.5)=nan;
index_tip  = index_tip(apple_p>0.5);
% apple    = movmean(apple,9);

% nframe   = 1:size(raw_apple,1);
nframe = 30*60;
index_tip = index_tip(1:nframe);
apple     = apple(1:nframe);
nframes   = 1: nframe;
index_tip = index_tip(1:nframe);
plot(index_tip,nframes,apple,nframes)
hold on
stem(Edge,nframe)
plot(apple,nframes,'color','red');
touch    = 0;time_fwd=[]; time_back=[];
trial_count = 0;
time_fwd=[]; time_back=[];
p=0;q=0;

for i = 9:nframes(end)-1
    if index_tip(i,1) <= Edge && index_tip(i+1,1) >= Edge&&index_tip(i+5,1)>Edge
        time_fwd =[time_fwd,i];
        trial_count = trial_count+1;
        hold on
        p=p+1;
        text(Edge,i,[num2str(p),'\rightarrow'],'Color','red','FontSize',15)
        j=i;
        while 1   
            j=j+1
            if j>=nframe(end)
                break
            end
            if index_tip(j,1) >= Edge && index_tip(j+1,1) <= Edge
               time_back = [time_back,j];
%                 if j ==time_back(q)
%                    j=nan;
%                 end
                q=q+1;
                text(Edge,j,[num2str(q),'\leftarrow'],'Color','green','FontSize',15)
                break
            end   
        end
    end
    
%     if TIP2(i,1) >= Edge && TIP2(i+1,1) <= Edge%&&TIP2(i+10,1)<Edge% && sum(isnan(apple_start(i-8:i))) > 0
%         time_back = [time_back,i];
%     end     
end
t=(time_back-time_fwd);

[location] = inthepicture(raw,M);
index_x  = [raw(location,M.base(1)),raw(location,M.MCP2(1)),raw(location,M.PIP2(1)),...
           raw(location,M.DIP2(1)),raw(location,M.tip2(1))];
index_y  = [raw(location,M.base(2)),raw(location,M.MCP2(2)),raw(location,M.PIP2(2)),...
           raw(location,M.DIP2(2)),raw(location,M.tip2(2))];       
        
        
thumb_x  = [raw(location,M.base(1)),raw(location,M.MCP1(1)),raw(location,M.PIP1(1)),raw(location,M.tip1(1))];
thumb_y  = [raw(location,M.base(2)),raw(location,M.MCP1(2)),raw(location,M.PIP1(2)),raw(location,M.tip1(2))];
%hand_y   = [raw(:,M.base(2)),raw(:,M.MCP1(2)),raw(:,M.MIP1(2)),raw(:,M.tip1(2))...
            %raw(:,M.MCP2(2)),raw(:,M.PIP2(2)),raw(:,M.DIP2(2)),raw(:,M.tip2(2))];
plot(index_x,index_y,'marker','+');
set(gca,'YDir','reverse')
hold on 
plot(thumb_x(:,1),thumb_y(:,1),'marker','diamond')

function [location_hand] = inthepicture(raw,M)  % find the lines with all the likelyhood more than 1
    likelyhood = [raw(:,M.base(3)),raw(:,M.MCP1(3)),raw(:,M.PIP1(3)),raw(:,M.tip1(3))...
              raw(:,M.MCP2(3)),raw(:,M.PIP2(3)),raw(:,M.DIP2(3)),raw(:,M.tip2(3))];
    likelyhood(likelyhood<0.9) = 0;
    mutiply_results      = likelyhood(:,1).*likelyhood(:,2).*likelyhood(:,3).*likelyhood(:,4)...
    .*likelyhood(:,5).*likelyhood(:,6).*likelyhood(:,7).*likelyhood(:,8);
    location_hand = find(mutiply_results~=0);
end 


% find the relative columns for each traking piont
function [column]  = findcolum(Txt,A,B,C)
        switch nargin  
            case 2
                 Aa        = strcmp(Txt(1:end),A);
                 column(1) = find(Aa==1);
            case 3
                Aa        = strcmp(Txt(1:end),A);       
                Ba        = strcmp(Txt(1:end),B);
                column(1) = find(Aa==1);
                column(2) = find(Ba==1);
            otherwise
                Aa        = strcmp(Txt(1:end),A);       
                Ba        = strcmp(Txt(1:end),B);
                Ca        = strcmp(Txt(1:end),C);
                column(1) = find(Aa==1);
                column(2) = find(Ba==1);
                column(3) = find(Ca==1);
        end          
end