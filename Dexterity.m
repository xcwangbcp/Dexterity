close all
clear 
monkey_name = {'2','11','25','35','43','44','60','70','76','132','133','137','159','187','195'};
%Step1  read the data into the workspace
% [F.filename_raw_hand,pathname] = uigetfile('*.csv','Pick a hand tracking csv file to load in');
% [F.raw_hand,F.Txt_hand,~]        = xlsread([pathname,F.filename_raw_hand]);
% F.filename_raw_apple           = uigetfile('*.csv','Pick an apple tracking csv file to load in ');
% [F.raw_apple,F.Txt_apple,~]      = xlsread([pathname,F.filename_raw_apple]);
% raw_apple = table2array(readtable(filename_raw_apple));
monkey_name={'2','11','43','44','60','132','133','137','159','187'};
monkey_name={'2'};
monkey_num = length(monkey_name );
p2mm       = 3; %1mm=3pixles
for j= 1:monkey_num
list        = ls([char(monkey_name(j)) '-*-hand.csv']);
length_file = size(list);
length_file = length_file(1); 
slitE_rate  = [];
wandE_rate  = [];
drop_rate   = [];
RT_valid    = [];
a2e_distance= [];
for i= 1:length_file
%     RT_valid = [];
    F.filename_raw_hand         = list(i,:);
    F.filename_raw_apple        = [F.filename_raw_hand(1:end-9) '-apple.csv'];
    F
    [F.raw_hand,F.Txt_hand,~]   = xlsread(F.filename_raw_hand);
    [F.raw_apple,F.Txt_apple,~] = xlsread(F.filename_raw_apple);
    R                           = singlefile(F);
%     RT_valid                    = mean([RT_valid; R.RT_valid(R.RT_valid<1000)]);

%     erro_rate(i)    = R.erro_rate;
    slitE_rate(i)   = R.slitE_rate;
    wandE_rate(i)   = R.wandE_rate;
%     grisE_rate(i)   = R.grisE_rate; 
    drop_rate(i)    = R.drop_rate;
    RT_valid(i)     = R.RT_vali_mean;
    a2e_distance(i) = R.a2e_distance/p2mm;
end

filename         = [F.filename_raw_hand(1:end-12) '.xlsx'];
% monkey           = ['#' filename(1:3)];
% m_d              = filename(4:9);
varNames         = {'ErrorType','Session'};
ErrorType        = {'slitErroRate';'wandErroRate';'dropErroRate';'fetchTime(Valid trials:ms)';'distance(Apple-Edge:mm)'};

T = table(ErrorType,[round(slitE_rate,2);round(wandE_rate,2);round(drop_rate,2);...
round(RT_valid,2);round(a2e_distance,2)],'VariableNames',varNames);
writetable(T,filename)
% the following write the data into an excel
% array    ={'Monkey','M/D','E_total','E_slithit','E_wand','E_grisp','RT','E_drop','D_apple2edge';...
%     monkey, m_d,round(erro_rate(1),2),round(slitE_rate(1),2),round(wandE_rate(1),2),...
%     round(grisE_rate(1),2),round(RT(1),2),round(drop_rate(1),2),round(a2e_distance(1),2);...
%     monkey, m_d,erro_rate_mean,slitE_rate_mean,wandE_rate_mean,grisE_rate_mean,RT_mean,...
%     drop_rate_mean,a2e_distance_mean};
% % data =[erro_rate_mean,slitE_rate_mean,wandE_rate_mean,grisE_rate_mean,RT_mean];
% xlswrite(filename,array);
end

function R = singlefile(F)

Hand   = read_csv(F.Txt_hand,F.raw_hand,F.raw_apple,F.filename_raw_hand,F.Txt_apple);
apple  = apple_trace(Hand);
action = search_action(Hand,apple);

% precision grisp error
% Step 9  find out the 4th type error, which is determined by the distance
% between the apple and the joints of index and thumb
% erroGrisp_trialID       = intersect(erroGrisp_trialID,invalid_id);
invalidID   = action.invalid;
dropID      = action.dropid;
sum_D_V     = unique([invalidID;dropID]);%;R.erroWandeID])
valid_action     = delete_errotrials(action.action,sum_D_V);
% take out the invalid trials and 
errorSlitHit = slitHitError(valid_action,Hand);
erroWander   = wanderError(valid_action,Hand);
erroGrisp    = precisGrispError(valid_action,Hand);

 R           = savewhat(apple,erroGrisp,errorSlitHit,erroWander,action);
end

% erroWander_trialID      = intersect(erroWander_trialID,invalid_id);
function R    = savewhat(apple,erroGrisp,errorSlitHit,erroWander,action)

R.applenum    = length(apple.apple_start);
% R.erroGrispID = erroGrisp.ID';
R.invalidID   = action.invalid;
R.dropID      = action.dropid;
R.sum_D_V     = unique([R.invalidID;R.dropID]);%;R.erroWandeID])
id_action     = delete_errotrials(action.action,R.sum_D_V);
R.RT_valid    = (id_action(:,7)-id_action(:,6))*1000/60;% take out the invalid and droptrials
R.RT_vali_mean= mean(R.RT_valid);   

R.errorSlitID = errorSlitHit.ID'; 
R.erroWandeID = erroWander.ID';
R.errorID     = unique([R.erroWandeID;R.errorSlitID]);%;;R.erroGrispID])
R.RT_erro     = (action.action(R.errorID ,7)-action.action(R.errorID ,6))*1000/60;
R.RT_all      = (action.action(:,7)-action.action(:,6))*1000/60;% include all the trials


% [locs,~]      = find(bsxfun(@eq,action.action(:,1),R.dropID));
id_act_correct= delete_errotrials(action.action,R.errorID);
R.correID     = id_act_correct(:,1);
R.RT_correct  = round((id_act_correct(:,7)-id_act_correct(:,6))*1000/60);% in ms

valid_trials  = R.applenum-length(R.sum_D_V);
R.erro_num    = length(R.errorID);
R.slitE_rate  = errorSlitHit.num/valid_trials ;
R.wandE_rate  = erroWander.num/valid_trials ;
R.grisE_rate  = erroGrisp.num/valid_trials ;
R.erro_rate   = R.erro_num/valid_trials ;
R.drop_rate   = length(R.dropID)/R.applenum ;
R.id_action   = action.action;
R.a2e_distance= action.distance;    
R.savefile    = erroGrisp.filename(1:end-5);
save([R.savefile '.mat'],'R');
end



function Hand = read_csv(Txt_hand,raw_hand,raw_apple,filename_raw_hand,Txt_apple)
threshold_hand  = 0.00001;
threshold_apple = 0.2;

Hand.tip1       = findcolum(Txt_hand,'tip1_x','tip1_y','tip1_p');
Hand.PIP1       = findcolum(Txt_hand,'PIP1_x','PIP1_y','PIP1_p');
Hand.MCP1       = findcolum(Txt_hand,'MCP1_x','MCP1_y','MCP1_p');
Hand.Base       = findcolum(Txt_hand,'base_x','base_y','base_p');

Hand.tip2       = findcolum(Txt_hand,'tip2_x','tip2_y','tip2_p');
Hand.DIP2       = findcolum(Txt_hand,'DIP2_x','DIP2_y','DIP2_p');
Hand.PIP2       = findcolum(Txt_hand,'PIP2_x','PIP2_y','PIP2_p');
Hand.MCP2       = findcolum(Txt_hand,'MCP2_x','MCP2_y','MCP2_p');   

[Hand.thumb_tip_x,Hand.thumb_tip_y] = readout(Hand.tip1,threshold_hand,raw_hand);
[Hand.thumb_pip_x,Hand.thumb_pip_y] = readout(Hand.PIP1,threshold_hand,raw_hand);
[Hand.thumb_mcp_x,Hand.thumb_mcp_y] = readout(Hand.MCP1,threshold_hand,raw_hand);
[Hand.base_x,          Hand.base_y] = readout(Hand.Base,threshold_hand,raw_hand);

[Hand.index_tip_x,Hand.index_tip_y] = readout(Hand.tip2,threshold_hand,raw_hand);
[Hand.index_dip_x,Hand.index_dip_y] = readout(Hand.DIP2,threshold_hand,raw_hand);
[Hand.index_pip_x,Hand.index_pip_y] = readout(Hand.PIP2,threshold_hand,raw_hand);
[Hand.index_mcp_x,Hand.index_mcp_y] = readout(Hand.MCP2,threshold_hand,raw_hand);

Hand.slot                   = findcolum(Txt_apple,'SlotTop_x','SlotBottom_x','SlotTop_y','SlotBottom_y');
Hand.edge_x                 = (mean(raw_apple(:,Hand.slot(1)))+mean(raw_apple(:,Hand.slot(2))))/2;
Hand.apple                  = findcolum(Txt_apple,'Apple_x','Apple_y','Apple_p');
[Hand.apple_x,Hand.apple_y] = readout(Hand.apple,threshold_apple,raw_apple);
Hand.filename_raw_hand      = filename_raw_hand;
M.pole                      = findcolum(Txt_apple,'Pole_x','Pole_y','Pole_p');
[Hand.pole_x,Hand.pole_y]   = readout(M.pole,threshold_apple,raw_apple);
% Hand.pole_y    = raw_apple(:,M.pole(2));
% Hand.pole_x    = raw_apple(:,M.pole(1));
end

function action = search_action(Hand,apple)
apple_start   = apple.apple_start;
% diff_apple    = apple.diff_apple;
apple_x       = apple.apple_x;
apple_y       = apple.apple_y;
pole_y        = mean(Hand.pole_y,'omitnan');
edge_x        = Hand.edge_x;
diff_polex    = diff(Hand.pole_x);
% diff_polex    = abs([nan;diff_polex]);
diff_apple    = diff(apple_x);
diff_apple    = [nan;diff_apple];
diff_apple    = abs(diff_apple);
apple_num     = length(apple_start);
drop_id       = [];
id_action     = zeros(apple_num,9);
% mostpositon   = apple.mostposition;
for j=1:apple_num
    if j<apple_num
        time_window = [apple_start(j):apple_start(j+1)-1];
    else 
        time_window = [apple_start(j):length(apple_x)];% the last apple show up and wait for 20 seconds
    end
    apple_window    = apple_x(time_window);
    pks = mean(findpeaks(apple_window));
    id_action (j,1) = j; % trial
    id_action (j,2) = apple_start(j); % apple come out time
%     diff_apple_window=diff(time_window);
    stop_loc       = continu_count_num(diff_apple(time_window),0.5,3);
    n=1;
    if length(stop_loc)>3
        apple_end = stop_loc(n);
        while apple_x(time_window(1)+ apple_end-1)>pks+20&&n<=length(stop_loc)-1 
           n=n+1;  
           apple_end = stop_loc(n);
        end
    else
        stop_loc  = find(diff_apple(time_window)<0.2,1,'first'); % -0.05
        if stop_loc
            apple_end = stop_loc;
        else
            apple_end = find(diff_apple(time_window)==min(abs(diff_apple(time_window))));
        end
    end 
   
    id_action (j,3) = time_window(1)+ apple_end-1; % apple stop time
    disappear       = continu_count_nan(apple_window,12); % 12 cotinuous Nans
%     disappear      = find(isnan(apple_window));% 
    m=1;
    if isempty(disappear)
        apple_disappear = apple_end+200;
    else
        apple_disappear = disappear(m);
        while apple_disappear<apple_end&&m<=length(disappear)-1
          m=m+1;
          apple_disappear = disappear(m);
        end
    end
%     if apple_disappear-stop_loc<20
%         apple_disappear = stop_loc+100;% sometimes the data make the disappear time close to the stoploc   
%     end
%     apple_disappear = find(diff_pole(time_window)>0.1,1,'first');% 
    if time_window(1)+ apple_disappear>length(apple_x)
       id_action (j,4) = length(apple_x)-10;
    else
       id_action (j,4)  = time_window(1)+ apple_disappear-2;% apple disappear time  
    end
    disappeartrace_y   = Hand.apple_y(id_action (j,2):time_window(end));
    locs               = find(disappeartrace_y>pole_y+6);%45 pixels=1.5cm
    diff_dropy         = diff(Hand.apple_y(locs));
    if length(locs)>2&&mean(diff_dropy)>0.0001
       drop_id         = [drop_id;j];
       drop_status     = 1;
    else
       drop_status     = 0;
    end
    id_action(j,9) = drop_status;
%     fprintf(fid,'%6d %6d %6d %6d %6d\n',apple_start_end_disappear(j,:));
end  
plot(apple_x(id_action(:,3)),id_action(:,3),'go');
plot(apple_x(id_action(:,4)),id_action(:,4),'b*')
hold off
stop_loc      = id_action(:,3);
stop_x_mean   = mean(apple_x(stop_loc));
stop_y_mean   = mean(apple_y(stop_loc));
window        = 30;% in pixle 1cm
for j=1:apple_num
    %if stop position is left or right 2 sigma than the mean,define it as invalid 
    m      = id_action(j,3);
    stop_x = mean(apple_x(m-2:m+2));
    stop_y = mean(apple_y(id_action(j,3)-2:id_action(j,3)+2));
    if stop_x>stop_x_mean+window ||...
       stop_x<stop_x_mean-window ||...
       stop_y>stop_y_mean+window ||...
       stop_y<stop_y_mean-window      
   
       invalid_status      = 1;
    else
       invalid_status      = 0;
    end
       id_action(j,5) = invalid_status;  
end

% invalid_id = unique(invalid_id);
% Step6:find out the errotypeI, which is grab the outside of the slitwant_count = 0; success_count = 0;
fwd_time       = []; 
back_time      = [];
% [id_action]    = delete_errotrials(id_action,invalid_id);
trials         = length(id_action(:,1));
for i= 1:trials
        for j=id_action(i,2):id_action(i,4)% from the show up time to apple_disappear_time(i)  
            if Hand.index_tip_x(j-1,1)<=Hand.edge_x...
               && Hand.index_tip_x(j+1,1)>=Hand.edge_x&&~isnan(apple_x(j))  %index_tip_x(j-2,1)<index_tip_x(j-1,1)&&&&index_tip_x(j+2,1)>index_tip_x(j+1,1)                 
               fwd_time  = [fwd_time,j];                              
               k = j;                                      
               while k<j+1000&&k<length(apple_x)                                 
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
        if isempty(fwd_time)||isempty(back_time)||back_time(end)<fwd_time(end)
            id_action(i,6:7) = [nan,nan];
            id_action(i,5)   = 1;%  sometimes the index never overpass the the edge
        else
            id_action(i,6:7) = [fwd_time(end); back_time(end)];% pass the inde
            index_Y          = Hand.index_tip_y(fwd_time(end):back_time(end));
            [~,touch_time ]  = max(index_Y );
            id_action(i,8)   = touch_time+fwd_time(end);
        end
        
end
action.action   = id_action;
action.invalid  = id_action((id_action(:,5)==1),1);
action.dropid   = drop_id ;
action.distance = stop_x_mean-edge_x;
end
function disappear_loc  = continu_count_nan(apple_window,Num)
         A              = apple_window;
         A(isnan(A))    = 1;
         B              = movsum(A,Num);
         disappear_loc  = find(B==Num);
end
% find out N=count number of apple coordinate smaller than the threshhold
function  stop_loc = continu_count_num(diff_apple_window,threshhold,count)
         A               = diff_apple_window;
         A(A<threshhold) = 1;
         B               = movsum(A,count);
         stop_loc        = find(B==count);         
end 
function apple = apple_trace(Hand)
apple_x = Hand.apple_x;
apple_y = Hand.apple_y;
% pole_y  = Hand.pole_y;
apple_x(apple_x<Hand.edge_x-20) = nan;
apple_x(apple_x>350)         = nan; % 340/350 is the board of the glass
% apple_x(apple_x>Hand.pole_x) = nan;
% apple_x(apple_y>Hand.pole_y) = nan;
% apple_y(isnan(apple_x))      = nan;
% apple_y(apple_y>pole_y+10)   = nan;
% apple_y(apple_y<pole_y-10)   = nan;
% apple_x(isnan(apple_y))      = nan;
backwards_time = 1*30;% backwards by human in 1 second
% Step 2:delet the part which apple is taken back by the pole by human 
diff_apple  = diff(apple_x);
diff_apple  = [nan;diff_apple];
sign_diff   = sign(diff_apple);% whether more or less than 0
sign_diff   = [nan;sign_diff];
sign_diff(isnan(sign_diff))=0;% change into 0 or +-1
sumwindow_sign   = movsum(sign_diff,backwards_time,'omitnan'); % 30 frames=1 seconds
locs_sign        = find(sumwindow_sign>=backwards_time);
[counts,centers] = hist(apple_x,20);
mostposition     = centers(counts==max(counts));
nframe           = length(Hand.index_tip_x);
nframes          = 1: nframe;
% delet the apple was taken back by human only the apple��s coordinates
% right than mostposition+20
if ~isempty(locs_sign )
    for i=1:length(locs_sign)-2
        if  apple_x(locs_sign(i))>mostposition+20
            apple_x(locs_sign(i)-2:locs_sign(i)+2)=nan;
        end 
    end
end
figure
plot(apple_x,nframes,'color','c')

%Step4: when the certain applecoming out,
%  find out when and the coordintes of the apple,
apple_start=[];% apple start means from the 1st appearance in the slit
for appleLength=3:length(apple_x)-5
    if isnan(apple_x(appleLength))&&isnan(apple_x(appleLength-1))&&isnan(apple_x(appleLength-2))
        if apple_x(appleLength+1)>apple_x(appleLength+2)&&apple_x(appleLength+2)>apple_x(appleLength+3)...
                &&apple_x(appleLength+3)>apple_x(appleLength+4)
            %if apple_x(appleLength+1)>280 % usually choose 280
                if apple_y(appleLength+1)>190&&apple_y(appleLength+1)<215 % in certain times, there are two apples were captured...
                    apple_start=[apple_start,appleLength+1];    %e.g. 76-7-7
                end
            %end
        end
    end
end
% the following used a more beutiful way to report when the apple start,but usually lost the 1st point.     
% point_start = find(z>30);
% diff_point_start = diff(point_start);
% % diff_point_start = [nan;diff_point_start];
% locs       = find(diff_point_start<10);%0.5s??????????????????��?????????10??
% point_start(locs)= nan;
% point_start=point_start(~isnan(point_start));
% for i=1:length(point_start)
%     if apple(point_start(i))<240
%         point_start(i)=nan;
%     end      
% end
hold on 
plot(apple_x(apple_start),apple_start,'ro');
hold on
plot(Hand.edge_x,nframes,'r')
title(Hand.filename_raw_hand(1:end-8));

% figure
% plot(pole_x,nframes,'color','c')
% fid=fopen([filename_raw_hand(1:end-4),'.txt'],'w');
% fprintf(fid,'%s %s %s %s %s\n','trialNo','appleStart','appleEnd','appleDisappear','grab/nograb');
% trial definition: a new/old apple come out 
apple.apple_start  = apple_start;
apple.diff_apple   = diff_apple;
apple.mostposition = mostposition;
apple.apple_x      = apple_x;
apple.apple_y      = apple_y;
end 

function erroGrisp   = precisGrispError(id_action,Hand)
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
 apple_x             = Hand.apple_x;
 apple_y             = Hand.apple_y;
 for i=1:trials
    time_window      = id_action(i,8)-2:id_action(i,8)+3;
    indexTip_window  = [index_tip_x(time_window),index_tip_y(time_window)];
    indexDip_window  = [index_dip_x(time_window),index_dip_y(time_window)];
    indexPip_window  = [index_pip_x(time_window),index_pip_y(time_window)];
    indexMcp_window  = [index_mcp_x(time_window),index_mcp_y(time_window)];
    
    thumbTip_window  = [thumb_tip_x(time_window),thumb_tip_y(time_window)];
    thumbPip_window  = [thumb_pip_x(time_window),thumb_pip_y(time_window)];
    thumbMcp_window  = [thumb_mcp_x(time_window),thumb_mcp_y(time_window)];
    appleP           = [apple_x(time_window),apple_y(time_window)];
    
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
%     title([filename_raw_hand(1:end-9)  ' trial #  ' num2str(id_action(i,1)) 'joints distance to apple(in pixle)'])
%     subplot(212)
%     plot(time_window,theta,'r-o')
%     title([filename_raw_hand(1:end-9)  ' trial #  ' num2str(id_action(i,1)) 'angle bet thumb and the index(in degree)'])
    if mean(d_indtip2thutip(:)>30)||mean(theta(:))>30
        erroGrisp_num     = erroGrisp_num+1;
        erroGrisp_trialID = [erroGrisp_trialID,id_action(i,1)];
    end
 end
 erroGrisp.num = erroGrisp_num;
 erroGrisp.ID  = erroGrisp_trialID;
 erroGrisp.filename = filename_raw_hand(1:end-4);
end

% Step 7:find out the error typeII,which is to pick out the wandering trials in the left
% If both fingers touched the food but the monkey released it and tried to pick it up again,
% we judged it as a ��wandering error��(nature 2012)
function erroWander = wanderError(id_action,Hand)
erroWander_num     = 0; 
erroWander_trialID = [];
trials             = length(id_action(:,1));

index_tip_x        = Hand.index_tip_x;
index_tip_y        = Hand.index_tip_y;
thumb_tip_x        = Hand.thumb_tip_x;
thumb_tip_y        = Hand.thumb_tip_y;
Edge_x             = Hand.edge_x;
apple_x            = Hand.apple_x;
filename_raw_hand  = Hand.filename_raw_hand;
index_tip_x(index_tip_x<Edge_x) = nan;
index_tip_x(index_tip_x>340)    = nan;
index_tip_y(isnan(index_tip_x)) = nan;
index_tip                       = [index_tip_x,406-index_tip_y];
thumb_tip                       = [thumb_tip_x,406-thumb_tip_y];

for j=1:trials
    time_window      = id_action(j,6):id_action(j,7);% which is the in+out of the slit
    apple_window     = apple_x(time_window);
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
erroWander.num = erroWander_num;
erroWander.ID  = erroWander_trialID;
end

% v_index= 
function errorSlitHit = slitHitError(id_action,Hand)
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
errorSlitHit.num = errorSlitHit_num;
errorSlitHit.ID  = errorSlitHit_trialID;
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

function id_action_validtrial = delete_errotrials(id_action,erroI_trialID)
        A=zeros(length(erroI_trialID),1);
        for i= 1:length(erroI_trialID)
            A(i)=find(id_action(:,1)==erroI_trialID(i));
        end
        id_action_validtrial=id_action;
        id_action_validtrial(A,:) = [];
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