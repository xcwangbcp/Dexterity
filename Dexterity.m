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
threshold  = 0.2;
index_tip(index_tip_p<threshold)=nan;
nframe    = length(index_tip);
nframes   = 1: nframe;
index_tip = movmean(index_tip,3);
% plot(index_tip,nframes,'color','red')


M.apple   = findcolum(Txt_apple,'Apple_x','Apple_y','Apple_p');
apple     = raw_apple(:,M.apple(1));
apple     = movmean(apple,5);% average in a second
apple_p   = raw_apple(:,M.apple(3));
% deletet the data by 3 criterier 
apple(apple_p<threshold)     = nan;
apple(apple<Edge )           = nan;
% arti_locs   = find(abs(diff_apple)>15);
% apple(arti_locs+1)=nan;
% delet the part which apple is taken back by the pole by human 
diff_apple  = diff(apple);
diff_apple  = [nan;diff_apple];
sign_diff   = sign(diff_apple);
sign_diff   = [nan;sign_diff];
sign_diff(isnan(sign_diff))=0;
sumwindow_sign = movsum(sign_diff,6,'omitnan');
locs_sign   = find(sumwindow_sign>=4);
for i=1:length(locs_sign)
    if apple(locs_sign(i))>240
       apple(locs_sign(i)-2:locs_sign(i)+2)=nan;
    end
end
plot(apple,nframes,'color','c')
apple_num = 0; 

x = apple(~isnan(apple));
p = ~isnan(apple); 
y = diff(x);
y =[nan;y];
z = NaN(1,nframe);
z(p)=y;
point = find(z>30); % 两帧之间的距离大于10个像素
diff_point = diff(point);
% diff_point = [nan;diff_point];
locs       = find(diff_point<10);%0.5s加一个约束条件，相邻的位置需要大于10帧
point(locs)= nan;
point=point(~isnan(point));
for i=1:length(point)
    if apple(point(i))<240
        point(i)=nan;
    end 
end
point=point(~isnan(point));
hold on
plot(apple(point),point,'marker','*','color','blue')
hold on
stem(Edge,nframe)

fwd_count = 0; want_count = 0; success_count = 0;
fwd_time  =[]; back_time  = [];
for i=1:length(point)-1
    time_window = [point(i),point(i+1)];
    hold on
    plot(index_tip(point(i):point(i+1)),point(i):point(i+1),'r')
    for j=point(i):point(i+1)                             % 从第i个苹果出来，到第i+1个苹果出来前的5帧内
       
         if index_tip(j,1)<= Edge && index_tip(j+1,1)>=Edge && sum(isnan(apple(j-8:j))) < 4  % 如果食指在该时间窗内伸入挡板，计为一次伸手
            fwd_count = fwd_count+1;
            fwd_time  = [fwd_time,j];                        %并且伸手的时候苹果在，即为想去拿苹果
%             if sum(isnan(apple(j-8:j)))>0
%                 want_count = want_count+1;
                k = j;                                        %进挡板后至少5帧，去判断出手的时间
                while 1
                    if index_tip(k,1)>=Edge&&index_tip(k+1,1)<=Edge&&sum(isnan(apple(k-8:k)))>0
                        success_count = success_count+1;
                        back_time     = [back_time,k];
                    end
                        k=k+1;
                    if  k>j+60
                        break
                    end
                end
%             end
         end
    end
end

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