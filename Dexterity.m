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
M.tip2       = findcolum(Txt_hand,'tip2_x','tip2_y','tip2_p');  
index_tip_x  = raw_hand(:,M.tip2(1));
index_tip_y  = raw_hand(:,M.tip2(2));
index_tip_p  = raw_hand(:,M.tip2(3));
threshold_hand  = 0.001;
index_tip_x(index_tip_p<threshold_hand)=nan;
index_tip_y(index_tip_p<threshold_hand)=nan;
index_tip_x = movmean(index_tip_x,5);
index_tip_y = movmean(index_tip_y,5);
% plot(index_tip,nframes,'color','red')
M.tip1       = findcolum(Txt_hand,'tip1_x','tip1_y','tip1_p');
thumb_tip_x  = raw_hand(:,M.tip1(1));
thumb_tip_y  = raw_hand(:,M.tip1(2));
thumb_tip_p  = raw_hand(:,M.tip1(3));
thumb_tip_x(thumb_tip_p<threshold_hand)=nan;
thumb_tip_x = movmean(thumb_tip_x,5);
thumb_tip_y(thumb_tip_p<threshold_hand)=nan;
thumb_tip_y = movmean(thumb_tip_y,5);

M.apple   = findcolum(Txt_apple,'Apple_x','Apple_y','Apple_p');
apple     = raw_apple(:,M.apple(1));
apple_y   = raw_apple(:,M.apple(2));
apple_p   = raw_apple(:,M.apple(3));
% M.pole    = findcolum(Txt_apple,'Pole_x','Plole_y');
% pole_y    = raw_apple(:,M.pole(2));
% deletet the data by 3 criterier 
threshold_apple = 0.1;
apple_y(apple_p<threshold_apple) = nan;
apple(apple_p<threshold_apple)   = nan;
apple(apple<Edge_x )             = nan;
apple(apple>340)                 = nan; % 350 is the board of the glass
apple     = movmean(apple,6);% average in 5 frames
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
nframe    = length(index_tip_x);
nframes   = 1: nframe;
% delet the apple was taken back by human only the apple‘s coordinates
% right than mostposition+20
% for i=1:length(locs_sign)-2
%     if apple(locs_sign(i))>mostposition+20
%        apple(locs_sign(i)-2:locs_sign(i)+2)=nan;
%     end
% end
plot(apple,nframes,'color','c')

%Step4: when the certain applecoming out,
%  find out when and the coordintes of the apple,
apple_start=[];% apple start means from the 1st appearance in the slit
for i=3:length(apple)-4
    if isnan(apple(i))&&isnan(apple(i-1))&&isnan(apple(i-2))
        if apple(i+1)>apple(i+2)&&apple(i+2)>apple(i+3)%&&apple(i+3)>apple(i+4)
            if apple(i+1)>260
                if apple_y(i+1)>195&&apple_y(i+1)<210 % in certain times, there are two apples were captured...
                    apple_start=[apple_start,i+1];    %e.g. 76-7-7
                end
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
hold off
nograb_num = 0;nograb_id = [];grab_status =[];
grab_id    =[];
fid=fopen([filename_raw_hand(1:end-4),'.txt'],'w');
fprintf(fid,'%s %s %s %s %s\n','trialNo','appleStart','appleEnd','appleDisappear','grab/nograb');
% Step 5: find out the no-grab trial and delete it ,make it related to the apple
% trial definition: a new/old apple come out 

apple_num = length(apple_start);
apple_start_end_disappear = zeros(apple_num,6);
id_action = zeros(apple_num,4);
for i=1:apple_num
    if i<apple_num
        time_window = [apple_start(i):apple_start(i+1)];
    else 
        time_window = [apple_start(i):length(apple)];% the last apple show up and wait for 20 seconds
    end
    apple_window     = apple(time_window);
    index_tip_x_window = index_tip_x(time_window);
    index_tip_y_window = index_tip_y(time_window);
%     figure    
%     figure    
%     plot(index_tip_x_window ,time_window,'r')
%     hold on 
%     plot(apple_window,time_window,'b-')
%     title(num2str(i));
% %     hold on
%     stem(Edge_x,nframe)
    apple_start_end_disappear(i,1) = i ; % trial
    apple_start_end_disappear(i,2) = apple_start(i); % apple come out time
    apple_end                      = find(diff_apple(time_window)>-0.02,1,'first');
    apple_start_end_disappear(i,3) = time_window(1)+ apple_end-1; % apple stop time
    apple_disappear                = find(isnan(apple_window),1,'first');% 判断这期间苹果在不在，不在的时间
    apple_start_end_disappear(i,4) = time_window(1)+ apple_disappear-2; % apple disappear time
    apple_start_end_disappear(i,6) = apple(apple_start_end_disappear(i,3));% 
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
id_action(nograb_id ,:) = [];
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
            if index_tip_x(j-1,1)<index_tip_x(j,1)&&index_tip_x(j,1)<= Edge_x && index_tip_x(j+1,1)>=Edge_x...
                 &&~isnan(apple(j))  %index_tip_x(j-2,1)<index_tip_x(j-1,1)&&&&index_tip_x(j+2,1)>index_tip_x(j+1,1)
             % 如果食指在该时间窗内伸入挡板，计为一次伸??,并且苹果??
                fwd_count = fwd_count+1;
                fwd_time  = [fwd_time,j];                   % 并且伸手的时候苹果在，即为想去拿苹果            
                
                k = j;                                      %进挡板后至少5帧，去判断出手的时间
                while 1                                      % 如果条件设成5点的，会损失某些??
                    if index_tip_x(k+1,1)<index_tip_x(k,1)&&index_tip_x(k,1)<=Edge_x...
                       &&index_tip_x(k-1,1)>=Edge_x       %&&index_tip_x(k-1,1)<index_tip_x(k-2,1) index_tip_x(k+2,1)<index_tip_x(k+1,1)&&
                        back_time     = [back_time,k];
                        j=k+1;
                        fprintf(fid,'%6d %6d %6d %6d \n',i,fwd_count,fwd_time(end),back_time(end));
                        break
                    end
                        k=k+1;
                end
            end
        end
        if  fwd_count==1
            correct_trialID = [correct_trialID,id_action(i,1)];
            correct_        = [fwd_time(end); back_time(end)];
%             correct_trial = [correct_trial',correct_]';
            id_action(i,6:7)=  correct_;
            index_Y         = index_tip_y(fwd_time(end):back_time(end));
            [~,touch_time ] = max(index_Y );
            id_action(i,8)  = touch_time+fwd_time(end);
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
fclose(fid);
% Step9,find out the 3rd type of error,which is hit the slit
erroIII_num =0; erroIII_trialID=[];
trials      = length(id_action(:,1));
for i=1:trials
    time_window  = id_action(i,6)-10:id_action(i,6)+10;
    v_index_x    = diff(index_tip_x(time_window));
    if v_index_x(10)<5
        erroIII_num     = erroIII_num+1;
        erroIII_trialID = [erroIII_trialID,id_action(i,1)];
    end
end
B=zeros(length(erroIII_trialID),1);
for i= 1:length(erroIII_trialID)
    B(i)=find(id_action(:,1)==erroIII_trialID(i));
end
id_action(B,:) = [];
% Step 7:find out the error typeII,which is to pick out the wandering trials in the left
% If both fingers touched the food but the monkey released it and tried to pick it up again,
% we judged it as a ‘wandering error’(nature 2012)
erroII_num =0; erroII_trialID=[];
trials      = length(id_action(:,1));
index_tip_x(index_tip_x<Edge_x)=nan;
index_tip_x(index_tip_x>340)   =nan;
index_tip_y(isnan(index_tip_x))=nan;
index_tip   = [index_tip_x,406-index_tip_y];
thumb_tip   = [thumb_tip_x,406-thumb_tip_y];
for i=1:trials
%     diff = id_action(i,3)-id_action(i,4);
%     if diff>1000
        time_window  = id_action(i,6):id_action(i,7);% which is the in+out of the slit
%     else
%         time_window  = id_action(i,2)-20:id_action(i,3);
%     end
    apple_window     = apple(time_window);
    index_tip_window = index_tip(time_window,:);
    thumb_tip_window = thumb_tip(time_window,:);
    figure 
    subplot(221)
    plot(index_tip_window(:,1),time_window,'r*')
    hold on 
    plot(apple_window,time_window,'b-')
    hold on
%     plot(index_tip_window(:,1),time_window,'r*')
    plot(thumb_tip_window(:,1),time_window,'g+')
    title([filename_raw_hand(1:end-8)  ' trial  ' num2str(id_action(i,1)) ' -x of index-thumb-apple'])
    hold off
    subplot(222)
    plot(index_tip_window(:,2),time_window,'r*')
    hold on 
    plot(thumb_tip_window(:,2),time_window,'g+')
    title([filename_raw_hand(1:end-8)  ' trial  ' num2str(id_action(i,1)) ' -y of index-thumb'])
%     set(gca,'YDir','reverse')
    hold off
    distance = sqrt((index_tip_window(:,1)-thumb_tip_window(:,1)).^2+(index_tip_window(:,2)-thumb_tip_window(:,2)).^2);
    subplot(223)
    plot(distance,time_window,'r')
    title([ ' trial  ' num2str(id_action(i,1)) ' Distance between index tip and thumb tip'])
%    [pks,locs]=findpeaks(distance,time_window,'MinPeakDistance',5);
   [tmax,vmax,tmin,vmin] = extrem_num(distance,time_window');
%    title([ ' trial  ' num2str(id_action(i,1)) ' Distance between index tip and thumb tip'])

  if ~isempty(tmin)&&length(tmax)>=2
      t     = [tmin;tmax];
      v     = [vmin;vmax];
      value = [t,v]; 
      value_sort = sortrows(value,1);
      value_sort_diff = diff(value_sort);
      pks  = length(value_sort_diff(:,1));
      N=0;
      for pks_num = 1:pks
          if value_sort_diff(pks_num,1)>5&&abs(value_sort_diff(pks_num,2))>9
             N=N+1;   
          end 
      end
      if N>=2
         erroII_num     = erroII_num+1;
         erroII_trialID = [erroII_trialID,id_action(i,1)];
      end
          
%       if vmax(1)-vmin(1)>10||vmax(2)-vmin(2)>10
%          erroII_num     = erroII_num+1;
%          erroII_trialID = [erroII_trialID,id_action(i,1)];
%       end
  end 
end
% Step 8  find out the 3rd type error, which is determined by the distance
% between the apple and the 



erro_num     = erroII_num+erroI_num+erroIII_num;
erro_rate    = erro_num/(apple_num-nograb_num);
delta_time    = mean((correct_trial(:,2)-correct_trial(:,1))*1000/60); % in ms unit 
% 
% v_index= 


function [tmax,vmax,tmin,vmin] = extrem_num(p,f)
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
% figure 
% plot them on a figure
subplot(224)
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