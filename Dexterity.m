close all
clear 
%Step1  read the data into the workspace
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
% Step 2:delet the part which apple is taken back by the pole by human 
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
% delet the apple was taken back by human only the apple‘s coordinates
% right than mostposition+20
for i=1:length(locs_sign)-2
    if apple(locs_sign(i))>mostposition+20
       apple(locs_sign(i)-2:locs_sign(i)+2)=nan;
    end
end
plot(apple,nframes,'color','c')

%Step4: when the certain applecoming out,
%  find out when and the coordintes of the apple,
apple_start=[];% apple start means from the 1st appearance in the slit
for i=3:length(apple)-4
    if isnan(apple(i))&&isnan(apple(i-1))&&isnan(apple(i-2))
        if apple(i+1)>apple(i+2)&&apple(i+2)>apple(i+3)%&&apple(i+3)>apple(i+4)
            if apple(i+1)>260
               apple_start=[apple_start,i+1];
            end
        end
    end
end
% the following used a more beutiful way to report when the apple start,but usually lost the 1st point.     
% point_start = find(z>30);
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
hold on 
plot(apple(apple_start),apple_start,'ro');
 
nograb_num = 0;nograb_id = [];grab_status =[];
grab_id    =[];
fid=fopen([filename_raw_hand(1:end-4),'.txt'],'w');
fprintf(fid,'%s %s %s %s %s\n','trialNo','appleStart','appleEnd','appleDisappear','grab/nograb');
% Step 5: find out the no-grab trial and delete it ,make it related to the apple
% trial definition: a new/old apple come out 

apple_num = length(apple_start);
apple_start_end_disappear = zeros(apple_num,5);
id_action = zeros(apple_num,4);
for i=1:apple_num
    if i<apple_num
        time_window = [apple_start(i):apple_start(i+1)];
    else 
        time_window = [apple_start(i):length(apple)];% the last apple show up and wait for 20 seconds
    end
    apple_window     = apple(time_window);
    index_tip_window = index_tip(time_window);
%     figure    
%     plot(index_tip_window ,time_window,'r')
%     hold on 
%     plot(apple_window,time_window,'b-')
%     title(num2str(i));
% %     hold on
%     stem(Edge_x,nframe)
    apple_start_end_disappear(i,1) = i ;
    apple_start_end_disappear(i,2) = apple_start(i);
    apple_end                      = find(diff_apple(time_window)>-0.6,1,'first');
    apple_start_end_disappear(i,3) = time_window(1)+ apple_end-1;
    apple_disappear                = find(isnan(apple_window),1,'first');% 判断这期间苹果在不在，不在的时间
    apple_start_end_disappear(i,4) = time_window(1)+ apple_disappear-2;
    % if disappear location is 5 ps right of the stable apple,means
    %the apple is taken back by experimenter,设定5,sometimes make mistakes
    if apple(apple_start_end_disappear(i,4))>mostposition+6          
         nograb_num          = nograb_num+1;
         nograb_id           = [nograb_id;i];
         id_action(i,:)      = [i,0,0,apple_start_end_disappear(i,4)];
         grab_status         = 0;
    else
        grab_id              = [grab_id,i];
        id_action (i,:)      = [i,time_window(1),time_window(end),apple_start_end_disappear(i,4)];
        grab_status          = 1;
    end
    apple_start_end_disappear(i,5) = grab_status;
    fprintf(fid,'%6d %6d %6d %6d %6d\n',apple_start_end_disappear(i,:));
end  
id_action(nograb_id,:) = [];
id_action(:,5:6) = 0;
% Step6:find out the errotypeI, which is grab the outside of the slitwant_count = 0; success_count = 0;
fwd_time  =[]; back_time  = [];erroI_trialID  = [];
correct_trial = [];correct_trialID = [];
erroI_num   = 0;
trials = length(grab_id);
for i= 1:trials
    fwd_count =0;
        for j=id_action(i,2)-20:id_action(i,3)%apple_disappear_time(i)  
            % 从第i个苹果出来，到第i个苹果消失的time window
            if index_tip(j-1,1)<index_tip(j,1)&&index_tip(j,1)<= Edge_x && index_tip(j+1,1)>=Edge_x...
                 &&~isnan(apple(j))  %index_tip(j-2,1)<index_tip(j-1,1)&&&&index_tip(j+2,1)>index_tip(j+1,1)
             % 如果食指在该时间窗内伸入挡板，计为一次伸??,并且苹果??
                fwd_count = fwd_count+1;
                fwd_time  = [fwd_time,j];                   % 并且伸手的时候苹果在，即为想去拿苹果            
                k = j;                                      %进挡板后至少5帧，去判断出手的时间
                while 1                                      % 如果条件设成5点的，会损失某些??
                    if index_tip(k+1,1)<index_tip(k,1)&&index_tip(k,1)<=Edge_x...
                       &&index_tip(k-1,1)>=Edge_x       %&&index_tip(k-1,1)<index_tip(k-2,1) index_tip(k+2,1)<index_tip(k+1,1)&&
                        back_time     = [back_time,k];
                        j=k+1;
                        fprintf(fid,'%6d %6d %6d %6d \n',i,fwd_count,fwd_time(end),back_time(end));
                        break
                    end
                        k=k+1;
                end
            end
        end
        if  fwd_count<=1
            correct_trialID = [correct_trialID,id_action(i,1)];
            correct_        = [fwd_time(end); back_time(end)];
%             correct_trial   = [correct_trial',correct_]';
            id_action(i,6:7)=  correct_;
        else
           erroI_num        = erroI_num+1;
           erroI_trialID     = [erroI_trialID,id_action(i,1)];
           id_action(i,5)   = 1;
           id_action(i,6:7)=  0;
        end
end
A=zeros(length(erroI_trialID),1);
for i= 1:length(erroI_trialID)
    A(i)=find(id_action(:,1)==erroI_trialID(i));
end
id_action(A,:) = [];
% Step 7:find out the error typeII,which is to pick out the wandering trials in the left
% If both fingers touched the food but the monkey released it and tried to pick it up again,
% we judged it as a ‘wandering error’(nature 2012)
erroII_num =0; erroII_trialID=[];
trials = length(id_action(:,1));
for i=1:trials
    time_window      = id_action(i,2)-20:id_action(i,3)-2*60;
    apple_window     = apple(time_window);
    index_tip_window = index_tip(time_window);
    index_tip_window(index_tip_window<Edge_x)=nan;
    index_tip_window(index_tip_window>260)=nan;
    apple_window     = apple_window(~isnan(index_tip_window));
    index_tip_window = index_tip_window(~isnan(index_tip_window));
    time_window      = time_window(~isnan(index_tip_window));
    figure    
    plot(index_tip_window ,time_window,'r')
    hold on 
    plot(apple_window,time_window,'b-')
    title(num2str(id_action(i,1)));
    
    [tmax,vmax,tmin,vmin] = extrem_num(index_tip_window,time_window');
end

erro_rate    = erro_num/(apple_num-nograb_num);
delta_time   = mean((correct_trial(:,2)-correct_trial(:,1))*1000/60); % in ms unit 
fclose(fid);
%grab_time = [fwd_time;back_time; back_time-fwd_time]';
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

function [tmax,vmax,tmin,vmin] = extrem_num(p,f);
v = p;
t = f;
Lmax = diff(sign(diff(v)))== -2; % logic vector for the local max value
Lmin = diff(sign(diff(v)))== 2; % logic vector for the local min value
% match the logic vector to the original vecor to have the same length
Lmax = [false; Lmax; false];
Lmin =  [false; Lmin; false];
tmax = t (Lmax); % locations of the local max elements
tmin = t (Lmin); % locations of the local min elements
vmax = v (Lmax); % values of the local max elements
vmin = v (Lmin); % values of the local min elements
 
% plot them on a figure
figure
plot(v,t);
xlabel('t'); ylabel('v');
hold on;
plot(vmax, tmax, 'r+');
plot(vmin,tmin, 'g+');
hold off;
end
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