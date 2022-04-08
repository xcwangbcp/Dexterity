close all
clear 
%  read the data into the workspace
[filename_raw_hand,pathname] = uigetfile('*.csv','Pick a hand tracking csv file to load in ');
[raw_hand,Txt_hand,~]        = xlsread([pathname,filename_raw_hand]);
filename_raw_apple           = uigetfile('*.csv','Pick an apple tracking csv file to load in ');
[raw_apple,Txt_apple,~]      = xlsread([pathname,filename_raw_apple]);
% filename_raw_hand  =  'D:\Code\Data\Apple\2-trimmed.csv';
% raw_apple = table2array(readtable(filename_raw_apple));
% 
M.slot     = findcolum(Txt_apple,'SlotTop_x','SlotBottom_x','SlotTop_y','SlotBottom_y');
Edge_x     = (mean(raw_apple(:,M.slot(1)))+mean(raw_apple(:,M.slot(2))))/2;
Edge_y_top = mean(raw_apple(:,M.slot(3)));
Edge_y_bot = mean(raw_apple(:,M.slot(4)));
% indexfinger tip
M.tip2     = findcolum(Txt_hand,'tip2_x','tip2_y','tip2_p');  
index_tip  = raw_hand(:,M.tip2(1));
index_tip_y= raw_hand(:,M.tip2(2));
index_tip_p= raw_hand(:,M.tip2(3));
threshold_hand  = 0.001;
index_tip(index_tip_p<threshold_hand)=nan;
index_tip = movmean(index_tip,5);

% plot(index_tip,nframes,'color','red')
M.apple   = findcolum(Txt_apple,'Apple_x','Apple_y','Apple_p');
apple     = raw_apple(:,M.apple(1));
% apple_y   = raw_apple(:,M.apple(2));
apple_p   = raw_apple(:,M.apple(3));
% deletet the data by 3 criterier 
threshold_apple = 0.2;
% apple_y(apple_p<threshold_apple)=nan;
apple(apple_p<threshold_apple)   = nan;
apple(apple<Edge_x )             = nan;
apple(apple>340)                 = nan; % 350 is the board of the glass
apple     = movmean(apple,5);% average in 5 frames
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
nframe    = length(index_tip);
nframes   = 1: nframe;
for i=1:length(locs_sign)-2
    if apple(locs_sign(i))>mostposition+20
       apple(locs_sign(i)-2:locs_sign(i)+2)=nan;
    end
end
plot(apple,nframes,'color','c')
apple_num = 0; 

locs_ss=[];
for i=3:length(apple)-4
    if isnan(apple(i))&&isnan(apple(i-1))&&isnan(apple(i-2))
        if apple(i+1)>apple(i+2)&&apple(i+2)>apple(i+3)%&&apple(i+3)>apple(i+4)
            if apple(i+1)>260
               locs_ss=[locs_ss,i+1];
            end
        end
    end
end
% the following used a more beutiful way to report when the apple start,but usually lost the 1st point.     
% point_start = find(z>30); % ????????????30??????
% diff_point_start = diff(point_start);
% % diff_point_start = [nan;diff_point_start];
% locs       = find(diff_point_start<10);%0.5s??????????????????λ?????????10??
% point_start(locs)= nan;
% point_start=point_start(~isnan(point_start));
% for i=1:length(point_start)
%     if apple(point_start(i))<240
%         point_start(i)=nan;
%     end 
% end
apple_start = locs_ss;% apple start means from the 1st appearance in the slit
hold on 
plot(apple(apple_start),apple_start,'ro');
fwd_count = 0; want_count = 0; success_count = 0;
fwd_time  =[]; back_time  = [];
apple_num = length(apple_start);
apple_start_end_disappear = zeros(apple_num,5); 
apple_start_end_disappear(:,1)=apple_start ;
nograb_num = 0;
fid=fopen([filename_raw_hand(1:end-4),'.txt'],'w');
for i=1:apple_num
    if i<apple_num
        time_window = [apple_start(i):apple_start(i+1)];
    else 
        time_window = [apple_start(i):length(apple)];% the last apple show up and wait for 20 seconds
    end
    apple_window     = apple(time_window);
    index_tip_window = index_tip(time_window);
    figure 
    title(num2str(i));
    plot(index_tip_window ,time_window,'r')
    hold on 
    plot(apple_window,time_window,'b-')
% %     hold on
%     stem(Edge_x,nframe)
    apple_end                      = find(diff_apple(time_window)>-0.02,1,'first');
    apple_start_end_disappear(i,2) = time_window(1)+ apple_end-1;
    apple_disappear                = find(isnan(apple_window),1,'first');% 判断这期间苹果在不在，不在的时间
    apple_start_end_disappear(i,3) = time_window(1)+ apple_disappear-2;% if disappear location is 10 ps right of the 
    apple_start_end_disappear(i,4) = i;
    apple_start_end_disappear(i,5) = nograb_num;
     fprintf(fid,'%6d %6d %6d %6d %6d\n',apple_start_end_disappear(i,:));
    if apple(apple_start_end_disappear(i,3))>mostposition+5           % stable apple,means, the apple is taken back by experimenter
         nograb_num=nograb_num+1;
    else
        for j=time_window(1):time_window(end)%apple_disappear_time(i)  % 从第i个苹果出来，到第i个苹果消失的time window
            if index_tip(j-2,1)<index_tip(j-1,1)&&index_tip(j-1,1)<index_tip(j,1)&&index_tip(j,1)<= Edge_x && index_tip(j+1,1)>=Edge_x...
                 &&index_tip(j+2,1)>index_tip(j+1,1)&&~isnan(apple(j))  
             % 如果食指在该时间窗内伸入挡板，计为一次伸??,并且苹果??
%              if index_tip_y(j)>130&&index_tip_y(j)<195    % if apple is no the pole,
                fwd_count = fwd_count+1;
                fwd_time  = [fwd_time,j];                   % 并且伸手的时候苹果在，即为想去拿苹果            
%             if sum(isnan(apple(j-8:j)))>0
%                 want_count = want_count+1;
                k = j;                                      %进挡板后至少5帧，去判断出手的时间
                while 1
                    if index_tip(k+1,1)<index_tip(k,1)&&index_tip(k,1)<=Edge_x...
                       &&index_tip(k-1,1)>=Edge_x       %&&index_tip(k-1,1)<index_tip(k-2,1) index_tip(k+2,1)<index_tip(k+1,1)&&
                        success_count = success_count+1;% 如果条件设成5点的，会损失某些??
                        back_time     = [back_time,k];
                        j=k+1;
                        fprintf(fid,' %6d %6d\n',fwd_time,back_time);
                        break
                    end
                        k=k+1;
                end
            end
        end   
    end
   
end
fclose(fid);
grab_time = [fwd_time;back_time; back_time-fwd_time]';
% delete the activity when the hand is out of the ROI which defined by the
hold on
plot(apple(apple_start_end_disappear(:,2)),apple_start_end_disappear(:,2),'go')
hold on
plot(apple(apple_start_end_disappear(:,3)-1),apple_start_end_disappear(:,3)-1,'b+')
j=1;
for  i=1:length(back_time)
    if index_tip_y(fwd_time(i))<130||index_tip_y(fwd_time(i))>200
        locs(j)=i;
        fwd_time(i)
        j=j+1;
    end 
end 
back_time(locs)=nan;fwd_time(locs)=nan;
grab_time = [ fwd_time(~isnan(fwd_time));back_time(~isnan(back_time))]';

% plot(thumb_x(:,1),thumb_y(:,1),'marker','diamond')
% find the lines with all the likelyhood more than 1
function [location_hand] = inthepicture(raw,M)  
    likelyhood = [raw(:,M.base(3)),raw(:,M.MCP1(3)),raw(:,M.PIP1(3)),raw(:,M.tip1(3))...
              raw(:,M.MCP2(3)),raw(:,M.PIP2(3)),raw(:,M.DIP2(3)),raw(:,M.tip2(3))];
    likelyhood(likelyhood<0.9) = 0;
    mutiply_results      = likelyhood(:,1).*likelyhood(:,2).*likelyhood(:,3).*likelyhood(:,4)...
    .*likelyhood(:,5).*likelyhood(:,6).*likelyhood(:,7).*likelyhood(:,8);
    location_hand = find(mutiply_results~=0);
end 


% find the relative columns for each traking piont
function [column]  = findcolum(Txt,A,B,C,D)
        switch nargin  
            case 2
                 Aa        = strcmp(Txt(1:end),A);
                 column(1) = find(Aa==1);
            case 3
                Aa        = strcmp(Txt(1:end),A);       
                Ba        = strcmp(Txt(1:end),B);
                column(1) = find(Aa==1);
                column(2) = find(Ba==1);
            case 4
                Aa        = strcmp(Txt(1:end),A);       
                Ba        = strcmp(Txt(1:end),B);
                Ca        = strcmp(Txt(1:end),C);
                column(1) = find(Aa==1);
                column(2) = find(Ba==1);
                column(3) = find(Ca==1);
            otherwise
                Aa        = strcmp(Txt(1:end),A);       
                Ba        = strcmp(Txt(1:end),B);
                Ca        = strcmp(Txt(1:end),C);
                Da        = strcmp(Txt(1:end),D);
                column(1) = find(Aa==1);
                column(2) = find(Ba==1);
                column(3) = find(Ca==1);
                column(4) = find(Da==1);
        end          
end