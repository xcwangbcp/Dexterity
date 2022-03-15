close all
clear 

[filename_raw_hand,pathname] = uigetfile('*.csv','Pick a hand tracking csv file to load in ');
[raw_hand,Txt_hand,~]        = xlsread([pathname,filename_raw_hand]);
filename_raw_apple           = uigetfile('*.csv','Pick an apple tracking csv file to load in ');
[raw_apple,Txt_apple,~]      = xlsread([pathname,filename_raw_apple]);
% filename_raw_hand  =  'D:\Code\Data\Apple\2-trimmed.csv';
% raw_apple = table2array(readtable(filename_raw_apple));
M.slot     = findcolum(Txt_apple,'SlotTop_x','SlotBottom_x');
Edge(:,1)  = mean(raw_apple(:,M.slot(1)));
Edge(:,2)  = mean(raw_apple(:,M.slot(2)));
Edge       = mean(Edge);
M.tip2     = findcolum(Txt_hand,'tip2_x','tip2_y','tip2_p'); % indexfinger tip 
index_tip  = raw_hand(:,M.tip2(1));
index_tip_p= raw_hand(:,M.tip2(3));
threshold_hand  = 0.001;
index_tip(index_tip_p<threshold_hand)=nan;
nframe    = length(index_tip);
nframes   = 1: nframe;
index_tip = movmean(index_tip,3);
% plot(index_tip,nframes,'color','red')


M.apple   = findcolum(Txt_apple,'Apple_x','Apple_y','Apple_p');
apple     = raw_apple(:,M.apple(1));
apple     = movmean(apple,5);% average in a second
apple_p   = raw_apple(:,M.apple(3));
% deletet the data by 3 criterier 
threshold_apple = 0.2;
apple(apple_p<threshold_apple) = nan;
apple(apple<Edge )             = nan;
% arti_locs   = find(abs(diff_apple)>15);
% apple(arti_locs+1)=nan;
% delet the part which apple is taken back by the pole by human 
diff_apple  = diff(apple);
diff_apple  = [nan;diff_apple];
sign_diff   = sign(diff_apple);% whether more or less than 0
sign_diff   = [nan;sign_diff];
sign_diff(isnan(sign_diff))=0;% change into 0 or +-1
sumwindow_sign   = movsum(sign_diff,6,'omitnan');
locs_sign        = find(sumwindow_sign>=4);
[counts,centers] = hist(apple,20);
mostposition = centers(counts==max(counts));
for i=1:length(locs_sign)
    if apple(locs_sign(i))>mostposition+5
       apple(locs_sign(i)-2:locs_sign(i)+2)=nan;
    end
end
plot(apple,nframes,'color','c')
apple_num = 0; 

x = apple(~isnan(apple));
p = ~isnan(apple); 
y = diff(x);
y =[nan;y];
% point_end = find(y>30)-10;
% if point_end(1)<=0
%     point_end(1) = 1;
% end
% m=0;
% while 1
%     m = m+1;
%     if isnan(y(point_end(m)))
%         y(point_end(m))=y(point_end(m+1));
%     else
%         break
%     end
% end

z = NaN(1,nframe);
% point_end=find(abs(z-y(point_end))<10^-4);
z(p)=y;
point_start = find(z>30); % 两帧之间的距离大于30个像素
diff_point_start = diff(point_start);
% diff_point_start = [nan;diff_point_start];
locs       = find(diff_point_start<10);%0.5s加一个约束条件，相邻的位置需要大于10帧
point_start(locs)= nan;
point_start=point_start(~isnan(point_start));
for i=1:length(point_start)
    if apple(point_start(i))<240
        point_start(i)=nan;
    end 
end
point_start=point_start(~isnan(point_start));
hold on
plot(apple(point_start),point_start,'marker','*','color','blue')


fwd_count = 0; want_count = 0; success_count = 0;
fwd_time  =[]; back_time  = [];
apple_num=length(point_start);
for i=1:apple_num
    if i<apple_num
        time_window = [point_start(i):point_start(i+1)];
    else 
        time_window = [point_start(i):length(apple)];% the last apple show up and wait for 20 seconds
    end
    apple_window     = apple(time_window);
    index_tip_window = index_tip(time_window);
    figure 
    plot(index_tip_window ,time_window,'r')
    hold on 
    plot(apple_window,time_window,'b-')
    hold on
    stem(Edge,nframe)
    apple_disappear         = find(isnan(apple_window),1,'first');% 判断这期间苹果在不在，不在的时间
    apple_disappear_time(i) = time_window(1)+ apple_disappear-1;

    for j=time_window(1):apple_disappear_time(i)                 % 从第i个苹果出来，到第i个苹果消失的time window
        
         if index_tip(j-2,1)<index_tip(j-1,1)&&index_tip(j-1,1)<index_tip(j,1)&&index_tip(j,1)<= Edge && index_tip(j+1,1)>=Edge...
                 &&index_tip(j+2,1)>index_tip(j+1,1)&&~isnan(apple(j))  % 如果食指在该时间窗内伸入挡板，计为一次伸手,并且苹果在
            fwd_count = fwd_count+1;
            fwd_time  = [fwd_time,j];                        %并且伸手的时候苹果在，即为想去拿苹果
%             if sum(isnan(apple(j-8:j)))>0
%                 want_count = want_count+1;
                k = j;                                        %进挡板后至少5帧，去判断出手的时间
                while 1
                    if index_tip(k+1,1)<index_tip(k,1)&&index_tip(k,1)<=Edge...
                       &&index_tip(k-1,1)>=Edge%&&index_tip(k-1,1)<index_tip(k-2,1) index_tip(k+2,1)<index_tip(k+1,1)&&
                        success_count = success_count+1;% 如果条件设成5点的，会损失某些点
                        back_time     = [back_time,k];
                        break
                    end
                        k=k+1;
%                     if  k>j+60
%                         break
%                     end
                end
%             end
         end
    end
end
grab_time = back_time-fwd_time;
locs      = grab_time>12;% 132-5-31 低于12帧的是没有获取到苹果的数据，其他猴子待验证
grab_time = grab_time(locs)/60;

for i = 9:nframes(end)-1
    if index_tip(i-2,1) < Edge &&index_tip(i-1,1) < Edge && index_tip(i,1) >= Edge&&index_tip(i+1,1)>Edge&&index_tip(i+2,1)>Edge
        if apple(i)>Edge&&apple(i)<230  % 20 pixles
        time_fwd =[time_fwd,i];
        trial_count = trial_count+1;
        hold on
        p=p+1;
        text(Edge-10,i,[num2str(p),'\rightarrow'],'Color','red','FontSize',15)
        j=i;
            while 1   
            j=j+1
                if j>=nframe(end)||j>=i+3*60;
                   break
                end
                if index_tip(j-2,1) > Edge&&index_tip(j-1,1) > Edge && index_tip(j,1) <= Edge&&index_tip(j+1,1)< Edge&&index_tip(j+2,1)<Edge
                    time_back = [time_back,j];
%                 if j ==time_back(q)
%                    j=nan;
%                 end
                    q=q+1;
                    text(Edge+10,j,[num2str(q),'\leftarrow'],'Color','green','FontSize',15)
                    break
                end   
            end
        end
    end
    
%     if TIP2(i,1) >= Edge && TIP2(i+1,1) <= Edge%&&TIP2(i+10,1)<Edge% && sum(isnan(apple_start(i-8:i))) > 0
%         time_back = [time_back,i];
%     end     
end









% for i= 5:length(apple)
%     if apple(i)>apple(i-1)&&apple(i+1)>apple(i)
%         
%     end
% end
hold on

% distance  = abs(diff(apple));
% apple_out(:,1) = apple(distance>30);
% apple_out(:,2) = find(distance>30)+2;
% % apple_out(:,1) = apple(apple_out(:,2));
% for i=1:length(apple_out(:,1))
%     if apple_out(i,1)>250
%         apple_out(i,1)=NaN;
%     end
% end
% plot(apple_out(:,1),apple_out(:,2),'marker','*','color','blue')

% pole
M.pole   = findcolum(Txt_apple,'Pole_x','Pole_y','Pole_p');
pole_x   = raw_apple(:,M.pole(1));
pole_p   = raw_apple(:,M.pole(3));
% nframe   = 1:size(raw_apple,1);
% nframe = 60*60;
index_tip = index_tip(1:nframe);
apple     = apple(1:nframe);
pole_x    = pole_x(1:nframe);
pole_x(pole_x<mode(pole_x) )     = nan;
pole_x(apple_p<threshold)= nan;
hold on
plot(pole_x,nframes,'c*')



 
% plot(pole_x,nframes,'color','blue')
touch    = 0; 
trial_count = 0;

p=0;q=0;


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