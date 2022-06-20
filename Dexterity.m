close all
clear 
%Step1  read the data into the workspace
[Hand.filename_raw_hand,pathname] = uigetfile('*.csv','Pick a hand tracking csv file to load in ');
[raw_hand,Txt_hand,~]        = xlsread([pathname,Hand.filename_raw_hand]);
filename_raw_apple           = uigetfile('*.csv','Pick an apple tracking csv file to load in ');
[raw_apple,Txt_apple,~]      = xlsread([pathname,filename_raw_apple]);
% filename_raw_hand  =  'D:\Code\Data\Apple\2-trimmed.csv';
% raw_apple = table2array(readtable(filename_raw_apple));
threshold_hand  = 0.00001;
threshold_apple = 0.1;

M.tip1          = findcolum(Txt_hand,'tip1_x','tip1_y','tip1_p');
M.PIP1          = findcolum(Txt_hand,'PIP1_x','PIP1_y','PIP1_p');
M.MCP1          = findcolum(Txt_hand,'MCP1_x','MCP1_y','MCP1_p');
M.Base          = findcolum(Txt_hand,'base_x','base_y','base_p');

M.tip2          = findcolum(Txt_hand,'tip2_x','tip2_y','tip2_p');
M.DIP2          = findcolum(Txt_hand,'DIP2_x','DIP2_y','DIP2_p');
M.PIP2          = findcolum(Txt_hand,'PIP2_x','PIP2_y','PIP2_p');
M.MCP2          = findcolum(Txt_hand,'MCP2_x','MCP2_y','MCP2_p');   

[Hand.thumb_tip_x,Hand.thumb_tip_y] = readout(M.tip1,threshold_hand,raw_hand);
[Hand.thumb_pip_x,Hand.thumb_pip_y] = readout(M.PIP1,threshold_hand,raw_hand);
[Hand.thumb_mcp_x,Hand.thumb_mcp_y] = readout(M.MCP1,threshold_hand,raw_hand);
[Hand.base_x,          Hand.base_y] = readout(M.Base,threshold_hand,raw_hand);

[Hand.index_tip_x,Hand.index_tip_y] = readout(M.tip2,threshold_hand,raw_hand);
[Hand.index_dip_x,Hand.index_dip_y] = readout(M.DIP2,threshold_hand,raw_hand);
[Hand.index_pip_x,Hand.index_pip_y] = readout(M.PIP2,threshold_hand,raw_hand);
[Hand.index_mcp_x,Hand.index_mcp_y] = readout(M.MCP2,threshold_hand,raw_hand);


M.slot      = findcolum(Txt_apple,'SlotTop_x','SlotBottom_x','SlotTop_y','SlotBottom_y');
Hand.edge_x = (mean(raw_apple(:,M.slot(1)))+mean(raw_apple(:,M.slot(2))))/2;
Edge_y_top  = mean(raw_apple(:,M.slot(3)));
Edge_y_bot  = mean(raw_apple(:,M.slot(4)));

M.apple           = findcolum(Txt_apple,'Apple_x','Apple_y','Apple_p');
[apple_x,apple_y] = readout(M.apple,threshold_apple,raw_apple);

% M.pole    = findcolum(Txt_apple,'Pole_x','Plole_y');
% pole_y    = raw_apple(:,M.pole(2));
% deletet the data by 3 criterier 
apple_x(apple_x<Hand.edge_x )        = nan;
apple_y(apple_x<Hand.edge_x )        = nan;
apple_x(apple_x>340)                 = nan; % 350 is the board of the glass
% Step 2:delet the part which apple is taken back by the pole by human 
diff_apple  = diff(apple_x);
diff_apple  = [nan;diff_apple];
sign_diff   = sign(diff_apple);% whether more or less than 0
sign_diff   = [nan;sign_diff];
sign_diff(isnan(sign_diff))=0;% change into 0 or +-1
sumwindow_sign   = movsum(sign_diff,6,'omitnan');
locs_sign        = find(sumwindow_sign>=3);
[counts,centers] = hist(apple_x,20);
mostposition     = centers(counts==max(counts));
nframe           = length(Hand.index_tip_x);
nframes          = 1: nframe;
% delet the apple was taken back by human only the appleï¿½ï¿½s coordinates
% right than mostposition+20
for i=1:length(locs_sign)-2
    if apple_x(locs_sign(i))>mostposition+20
       apple_x(locs_sign(i)-2:locs_sign(i)+2)=nan;
    end 
end
plot(apple_x,nframes,'color','c')

%Step4: when the certain applecoming out,
%  find out when and the coordintes of the apple,
apple_start=[];% apple start means from the 1st appearance in the slit
for appleLength=3:length(apple_x)-5
    if isnan(apple_x(appleLength))&&isnan(apple_x(appleLength-1))&&isnan(apple_x(appleLength-2))
        if apple_x(appleLength+1)>apple_x(appleLength+2)&&apple_x(appleLength+2)>apple_x(appleLength+3)&&apple_x(appleLength+3)>apple_x(appleLength+4)
            if apple_x(appleLength+1)>280
                if apple_y(appleLength+1)>190&&apple_y(appleLength+1)<215 % in certain times, there are two apples were captured...
                    apple_start=[apple_start,appleLength+1];    %e.g. 76-7-7
                end
            end
        end
    end
end
% the following used a more beutiful way to report when the apple start,but usually lost the 1st point.     
% point_start = find(z>30);
% diff_point_start = diff(point_start);
% % diff_point_start = [nan;diff_point_start];
% locs       = find(diff_point_start<10);%0.5s??????????????????ï¿½ï¿½?????????10??
% point_start(locs)= nan;
% point_start=point_start(~isnan(point_start));
% for i=1:length(point_start)
%     if apple(point_start(i))<240
%         point_start(i)=nan;
%     end      
% end
hold on 
plot(apple_x(apple_start),apple_start,'ro');
title(Hand.filename_raw_hand(1:end-8));
hold off
invalid_id = [];grab_status =[];
grab_id    =[];
% fid=fopen([filename_raw_hand(1:end-4),'.txt'],'w');
% fprintf(fid,'%s %s %s %s %s\n','trialNo','appleStart','appleEnd','appleDisappear','grab/nograb');
% Step 5: find out the no-grab trial and delete it ,make it related to the apple
% trial definition: a new/old apple come out 

apple_num                 = length(apple_start);
apple_start_end_disappear = zeros(apple_num,6);
id_action                 = zeros(apple_num,4);
for j=1:apple_num
    if j<apple_num
        time_window = [apple_start(j):apple_start(j+1)];
    else 
        time_window = [apple_start(j):length(apple_x)];% the last apple show up and wait for 20 seconds
    end
    apple_window       = apple_x(time_window);
    index_tip_x_window = Hand.index_tip_x(time_window);
    index_tip_y_window = Hand.index_tip_y(time_window);
%     figure    
%     plot(index_tip_x_window ,time_window,'r')
%     hold on 
%     plot(apple_window,time_window,'b-')
%     title(num2str(i));
    id_action (j,1) = j ; % trial
    id_action (j,2) = apple_start(j); % apple come out time
    apple_end       = find(diff_apple(time_window)>-0.05,1,'first');
    id_action (j,3) = time_window(1)+ apple_end-1; % apple stop time
    apple_disappear = find(isnan(apple_window),1,'first');% 
    id_action (j,4) = time_window(1)+ apple_disappear-2; % apple disappear time
%     apple_start_end_disappear(j,6) = apple_x(apple_start_end_disappear(j,3));% 
    % if disappear location is 5 ps right of the stable apple,means
    %the apple is taken back by experimenter,ï¿½è¶¨5,sometimes make mistakes
    if apple_x(id_action(j,4))>mostposition+8||...
       apple_y(id_action(j,4))>apple_y(id_action(j,3))+6
%        apple_y(id_action(j,4))<apple_y(id_action(j,3))-6          
       invalid_id          = [invalid_id;j];
       grab_status         = 0;
    else
        grab_id            = [grab_id,j];
        grab_status        = 1;
    end
    id_action(j,5) = grab_status;
%     fprintf(fid,'%6d %6d %6d %6d %6d\n',apple_start_end_disappear(j,:));
end  
invalid_id = unique(invalid_id);
% id_action(invalid_id ,:) = [];
% Step6:find out the errotypeI, which is grab the outside of the slitwant_count = 0; success_count = 0;
fwd_time       = []; 
back_time      = [];
erroI_trialID  = [];
erroI_num      = 0;
trials         = length(id_action(:,1));
for i= 1:trials
    fwd_count =0; % 
        for j=id_action(i,2):id_action(i,4)% from the show up time to apple_disappear_time(i)  
            if Hand.index_tip_x(j-1,1)<Hand.index_tip_x(j,1)&&Hand.index_tip_x(j,1)<=Hand.edge_x...
               && Hand.index_tip_x(j+1,1)>=Hand.edge_x&&~isnan(apple_x(j))  %index_tip_x(j-2,1)<index_tip_x(j-1,1)&&&&index_tip_x(j+2,1)>index_tip_x(j+1,1)
               fwd_count = fwd_count+1;                    % 
               fwd_time  = [fwd_time,j];                   %           
               k = j;                                      %
               while k<j+1000                                     % 
%                     if Hand.index_tip_x(k+1,1)<Hand.index_tip_x(k,1)&&Hand.index_tip_x(k,1)<=Hand.edge_x...
%                        &&Hand.index_tip_x(k-1,1)>=Hand.edge_x    %&&index_tip_x(k-1,1)<index_tip_x(k-2,1) index_tip_x(k+2,1)<index_tip_x(k+1,1)&&
                        if Hand.index_tip_x(k,1)>=Hand.edge_x&&Hand.index_tip_x(k+1,1)<=Hand.edge_x
                           back_time     = [back_time,k];
                           j=k+1;
%                         fprintf(fid,'%6d %6d %6d %6d \n',j,fwd_count,fwd_time(end),back_time(end));
                        break
                        end
                        k=k+1;
                end
            end
        end
         id_action(i,6:7) = [fwd_time(end); back_time(end)];% pass the inde
         index_Y          = Hand.index_tip_y(fwd_time(end):back_time(end));
         [~,touch_time ]  = max(index_Y );
         id_action(i,8)   = touch_time+fwd_time(end);
%         if  fwd_count==1
%             correct_trialID = [correct_trialID,id_action(i,1)];
%             correct_        = [fwd_time(end); back_time(end)];
% %             correct_trial = [correct_trial',correct_]';
%             id_action(i,6:7)=  correct_;
%             index_Y         = Hand.index_tip_y(fwd_time(end):back_time(end));
%             [~,touch_time ] = max(index_Y );
%             id_action(i,8)  = touch_time+fwd_time(end);
%         else
%            erroI_num        = erroI_num+1;
%            erroI_trialID     = [erroI_trialID,id_action(i,1)];
%            id_action(i,5)   = 1;
%            id_action(i,6:7)=  0;
%         end
end
% [id_action]=delete_errotrials(id_action,erroI_trialID);
% fclose(fid);

% Step9,find out the 3rd type of error,which is hit the slit
[erroGrisp_num,erroGrisp_trialID]   = precisGrispError(id_action,Hand,apple_x,apple_y);
% erroGrisp_trialID       = intersect(erroGrisp_trialID,invalid_id);
% erroGrisp_num           = length(erroGrisp_trialID);

[errorSlitHit_num,errorSlitHit_trialID]= slitHitError(id_action,Hand);
% errorSlitHit_trialID    = intersect(errorSlitHit_trialID,invalid_id);
% errorSlitHit_num        = length(errorSlitHit_trialID );

[erroWander_num ,erroWander_trialID] = wanderError(id_action,Hand,apple_x);
% erroWander_trialID      = intersect(erroWander_trialID,invalid_id);
% erroWander_num          = length(erroWander_trialID);

R.applenum    = apple_num;
R.erroGrispID = erroGrisp_trialID';
R.errorSlitID = errorSlitHit_trialID';
R.erroWandeID = erroWander_trialID';
R.invalidID   = invalid_id;
R.errorID     = unique ([errorSlitHit_trialID erroGrisp_trialID erroWander_trialID]);
R.RT_all      = round((id_action(:,7)-id_action(:,6))*1000/60);
[locs,~]      = find(bsxfun(@eq,id_action(:,1),R.errorID));
R.RT_error    = round((id_action(locs,7)-id_action(locs,6))*1000/60);
id_act_correct= delete_errotrials(id_action,R.errorID);
R.correID     = id_act_correct(:,1);
R.RT_correct  = round((id_act_correct(:,7)-id_act_correct(:,6))*1000/60);
% in ms
invalid_num    = length(invalid_id);
R.erro_num     = length(R.errorID);
R.slitE_rate   = errorSlitHit_num/apple_num;
R.wandE_rate   = erroWander_num/apple_num;
R.grisE_rate   = erroGrisp_num/apple_num;
R.erro_rate    = (R.erro_num-invalid_num)/(apple_num-invalid_num);
R.savefile     = [Hand.filename_raw_hand(1:end-9) '.mat'];
R.id_action    = id_action;
save(R.savefile,'R');
clear
% delta_time   = mean((correct_trial(:,2)-correct_trial(:,1))*1000/60); % in ms unit 
% precision grisp error
% Step 9  find out the 4th type error, which is determined by the distance
% between the apple and the joints of index and thumb
function [erroGrisp_num,erroGrisp_trialID] = precisGrispError(id_action,Hand,apple,apple_y)
 erroGrisp_num       = 0; 
 erroGrisp_trialID   = [];
 trials              = length(id_action(:,1));
 index_tip_x         = Hand.index_tip_x;
 index_tip_y         = Hand.index_tip_y;
 index_dip_x         = Hand.index_dip_x;
 index_dip_y         = Hand.index_dip_y;
 index_pip_x         = Hand.index_pip_x;
 index_pip_y         = Hand.index_pip_y;
 index_mcp_x         = Hand.index_mcp_x;
 index_mcp_y         = Hand.index_mcp_y;
 thumb_tip_x         = Hand.thumb_tip_x;
 thumb_tip_y         = Hand.thumb_tip_y;
 thumb_pip_x         = Hand.thumb_pip_x;
 thumb_pip_y         = Hand.thumb_pip_y;
 thumb_mcp_x         = Hand.thumb_mcp_x;
 thumb_mcp_y         = Hand.thumb_mcp_y;
 base_x              = Hand.base_x;
 base_y              = Hand.base_y;
 filename_raw_hand   = Hand.filename_raw_hand;
 for i=1:trials
    time_window      = id_action(i,8):id_action(i,8)+3;
    indexTip_window  = [index_tip_x(time_window),index_tip_y(time_window)];
    indexDip_window  = [index_dip_x(time_window),index_dip_y(time_window)];
    indexPip_window  = [index_pip_x(time_window),index_pip_y(time_window)];
    indexMcp_window  = [index_mcp_x(time_window),index_mcp_y(time_window)];
    
    thumbTip_window  = [thumb_tip_x(time_window),thumb_tip_y(time_window)];
    thumbPip_window  = [thumb_pip_x(time_window),thumb_pip_y(time_window)];
    thumbMcp_window  = [thumb_mcp_x(time_window),thumb_mcp_y(time_window)];
    appleP           = [apple(time_window),apple_y(time_window)];
    
    base_window      = [base_x(time_window), base_y(time_window)];
    L                = length(time_window);
    d_indextip2apple = zeros(L,1);
    d_indexpip2apple = zeros(L,1);
    d_indexdip2apple = zeros(L,1);
    d_thumbtip2apple = zeros(L,1);
    d_thumbpip2apple = zeros(L,1);
    d_indtip2thutip  = zeros(L,1);
    theta            = zeros(L,1);
    for j=1:L
        d_indextip2apple(j,1)  = norm(indexTip_window(j,:)-appleP(j,:));
        d_indexdip2apple(j,1)  = norm(indexDip_window(j,:)-appleP(j,:));
        d_indexpip2apple(j,1)  = norm(indexPip_window(j,:)-appleP(j,:));
        
        d_indtip2thutip(j,1)   = norm(indexTip_window(j,:)-thumbTip_window(j,:));
        
        d_thumbtip2apple(j,1)  = norm(thumbTip_window(j,:)-appleP(j,:));
        d_thumbpip2apple(j,1)  = norm(thumbPip_window(j,:)-appleP(j,:));
        aa = norm(base_window(j,:)-thumbMcp_window(j,:));
        bb = norm(thumbMcp_window(j,:)-indexMcp_window(j,:));
        cc = norm(base_window(j,:)-indexMcp_window(j,:));
        theta(j)= acosd((aa^2+cc^2-bb^2)/2/aa/cc);                   
    end
%     figure
%     subplot(211)
%     plot(time_window,d_indextip2apple,'r-o',time_window,d_indexdip2apple,'g-*',time_window,d_indexpip2apple,...
%         'b-s',time_window,d_thumbtip2apple,'k-^',time_window,d_thumbpip2apple,'y-v',time_window,d_indtip2thutip,...
%         'm-+')
%     legend('indextip2apple','indexdip2apple','indexpip2apple','thumbtip2apple','thumbpip2apple','indtip2thutip')
%     subplot(212)
%     plot(time_window,theta,'r-o')
%     title([filename_raw_hand(1:end-4)  ' trial #  ' num2str(id_action(i,1)) ' degree of »¢¿Ú(in pixle)'])
    if d_indtip2thutip(:)>30|theta(:)>30
        erroGrisp_num     = erroGrisp_num+1;
        erroGrisp_trialID = [erroGrisp_trialID,id_action(i,1)];
    end
 end
end



% Step 7:find out the error typeII,which is to pick out the wandering trials in the left
% If both fingers touched the food but the monkey released it and tried to pick it up again,
% we judged it as a ï¿½ï¿½wandering errorï¿½ï¿½(nature 2012)
function [erroWander_num ,erroWander_trialID] = wanderError(id_action,Hand,apple)
erroWander_num     = 0; 
erroWander_trialID = [];
trials             = length(id_action(:,1));

index_tip_x        = Hand.index_tip_x;
index_tip_y        = Hand.index_tip_y;
thumb_tip_x        = Hand.thumb_tip_x;
thumb_tip_y        = Hand.thumb_tip_y;
Edge_x             = Hand.edge_x;
filename_raw_hand  = Hand.filename_raw_hand;
index_tip_x(index_tip_x<Edge_x) = nan;
index_tip_x(index_tip_x>340)    = nan;
index_tip_y(isnan(index_tip_x)) = nan;
index_tip                       = [index_tip_x,406-index_tip_y];
thumb_tip                       = [thumb_tip_x,406-thumb_tip_y];

for j=1:trials
    time_window      = id_action(j,6):id_action(j,7);% which is the in+out of the slit
    apple_window     = apple(time_window);
    index_tip_window = index_tip(time_window,:);
    thumb_tip_window = thumb_tip(time_window,:);
    distance         = sqrt((index_tip_window(:,1)-thumb_tip_window(:,1)).^2+(index_tip_window(:,2)-thumb_tip_window(:,2)).^2);
%     figure
%     subplot(221)
%     plot(index_tip_window(:,1),time_window,'r*')
%     hold on 
%     plot(apple_window,time_window,'b-')
%     hold on
%     plot(thumb_tip_window(:,1),time_window,'g+')
%     title([filename_raw_hand(1:end-8)  ' trial  ' num2str(id_action(j,1)) ' -x of index-thumb-apple'])
%     hold off
%     subplot(222)
%     plot(index_tip_window(:,2),time_window,'r*')
%     hold on 
%     plot(thumb_tip_window(:,2),time_window,'g+')
%     title([filename_raw_hand(1:end-8)  ' trial  ' num2str(id_action(j,1)) ' -y of index-thumb'])
% %     set(gca,'YDir','reverse')
%     hold off
%     subplot(223)
%     plot(distance,time_window,'r')
%     title([ ' trial  ' num2str(id_action(j,1)) ' Distance between index tip and thumb tip'])
%    [pks,locs]=findpeaks(distance,time_window,'MinPeakDistance',5);
   [tmax,vmax,tmin,vmin] = extrem_num(distance,time_window');
  if ~isempty(tmin)&&length(tmax)>=2
      t     = [tmin;tmax];
      v     = [vmin;vmax];
      value = [t,v]; 
      value_sort      = sortrows(value,1);
      value_sort_diff = diff(value_sort);
      pks             = length(value_sort_diff(:,1));
      N               = 0;
      for pks_num = 1:pks
          if value_sort_diff(pks_num,1)>5&&abs(value_sort_diff(pks_num,2))>9
             N=N+1;   
          end 
      end
      if N>=2
         erroWander_num     = erroWander_num+1;
         erroWander_trialID = [erroWander_trialID,id_action(j,1)];
      end
%       if vmax(1)-vmin(1)>10||vmax(2)-vmin(2)>10
%          erroII_num     = erroII_num+1;
%          erroII_trialID = [erroII_trialID,id_action(i,1)];
%       end
  end 
end
end

% v_index= 
function [errorSlitHit_num,errorSlitHit_trialID] = slitHitError(id_action,Hand)
errorSlitHit_num     = 0;
errorSlitHit_trialID = [];
trials               = length(id_action(:,1));
filename_raw_hand    = Hand.filename_raw_hand;
index_tip_x          = Hand.index_tip_x;
Edge_x               = Hand.edge_x;
for i=1:trials
    time_window  = id_action(i,6)-10:id_action(i,6)+10;
    v_index_x    = [false;diff(index_tip_x(time_window))];
%     figure
%     subplot(211)
%     plot(index_tip_x(time_window),time_window,'r-*')
%     hold on 
%     plot([Edge_x,Edge_x],[time_window(1),time_window(end)],'g-')
%     plot([min(index_tip_x(time_window)),max(index_tip_x(time_window))],[id_action(i,6),id_action(i,6)],'g-')
%     title([filename_raw_hand(1:end-4)  ' trial #  ' num2str(id_action(i,1)) ' -x of index'])
%     hold off
%     subplot(212)
%     plot(v_index_x,time_window,'r-*')
%     hold on 
%     plot([min(v_index_x),max(v_index_x)],[id_action(i,6),id_action(i,6)],'g-')
%     title([filename_raw_hand(1:end-4)  ' trial #  ' num2str(id_action(i,1)) ' -x of index speed(in pixle/frame)'])
%     hold off
    if mean(v_index_x(9:11))<5 % 
        errorSlitHit_num     = errorSlitHit_num+1;
        errorSlitHit_trialID = [errorSlitHit_trialID,id_action(i,1)];
    end
end
% [id_action]    = delete_errotrials(id_action,errorSlitHit_trialID);
end 



function [x,y] = readout(M,threshold,raw_hand)
    x  = raw_hand(:,M(1));
    y  = raw_hand(:,M(2));
    p  = raw_hand(:,M(3));
    x(p<threshold)=nan;
    y(p<threshold)=nan;
    x = movmean(x,5);
    y = movmean(y,5);
end

function [id_action]=delete_errotrials(id_action,erroI_trialID)
        A=zeros(length(erroI_trialID),1);
        for i= 1:length(erroI_trialID)
            A(i)=find(id_action(:,1)==erroI_trialID(i));
        end
        id_action(A,:) = [];
end

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
%     subplot(224)
%     plot(v,t);
%     xlabel('t'); ylabel('v');
%     hold on;
%     plot(vmax, tmax, 'r+');
%     plot(vmin,tmin, 'g+');
%     hold off;
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