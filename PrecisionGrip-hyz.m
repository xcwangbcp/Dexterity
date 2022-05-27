function [precision_grip_error,error_rate3] = PrecisionGrip(fps,raw_hand,trial_count,time_TIP2,trial_start,apple_position)
TIP1 = movmean(raw_hand(:,10:11),9);
PIP1 = movmean(raw_hand(:,7:8),9);
MCP1 = movmean(raw_hand(:,4:5),9);
TIP2 = movmean(raw_hand(:,22:23),9);
DIP2 = movmean(raw_hand(:,19:20),9);
PIP2 = movmean(raw_hand(:,16:17),9);
MCP2 = movmean(raw_hand(:,13:14),9);
base = movmean(raw_hand(:,1:2),9);

lengthT1P1 = ones(size(TIP1,1),3)*nan;
distanceT1A = ones(size(TIP1,1),3)*nan;
th_MBM = ones(size(TIP1,1),3)*nan;

for i = 1:length(trial_start)
    for j = -3:4+time_TIP2(i)*fps
        k = trial_start(i)+j-1;
        lengthT1P1(k,:) = [norm(TIP1(k,:)-PIP1(k,:)),i,j];
        distanceT1A(k,:) = [norm(TIP1(k,:)-apple_position(i,:)),i,j];
        th_MBM(k,:) = [180/pi*acos(dot(MCP1(k,:)-base(k,:),MCP2(k,:)-base(k,:))/(norm(MCP1(k,:)-base(k,:))*norm(MCP2(k,:)-base(k,:)))),i,j];
    end    
end
lengthT1P1 = [movmean(lengthT1P1(:,1),9),lengthT1P1(:,2:3)];

th_MBM = [movmean(th_MBM(:,1),9),th_MBM(:,2:3)];
plot(th_MBM(:,1),'r')

th_min = [];
for i = 1:trial_count
    th_min = [th_min,min(th_MBM(find(th_MBM(:,2) == i),1))];
end

distanceT1A = [movmean(distanceT1A(:,1),9),distanceT1A(:,2:3)];
plot(distanceT1A(:,1),'b')
[~,trough_locs] = findpeaks(-distanceT1A(:,1),'MinPeakProminence',1);

team = [];
for i = 1:length(trough_locs)
    team = [team;find(trial_start < trough_locs(i),1,'last')];
end
trough_locs = [team,trough_locs];

precision_grip_error = [];
for i = 1:trial_count
    if isempty(find(trough_locs(:,1) == i)) == 0
        [trough_row,~] = find(trough_locs(:,1) == i);
        j = trough_locs(trough_row(end),2);
        thumb_apple = cross([TIP1(j,:)-PIP1(j,:),0],[apple_position(i,:)-PIP1(j,:),0]);
        indexfinger_apple = (abs(det([TIP2(j,:)-DIP2(j,:);apple_position(i,:)-DIP2(j,:)]))/norm(TIP2(j,:)-DIP2(j,:)))-(abs(det([PIP2(j,:)-DIP2(j,:);apple_position(i,:)-DIP2(j,:)]))/norm(PIP2(j,:)-DIP2(j,:)));
        if th_min(i) < 10
            precision_grip_error = [precision_grip_error,i];
        elseif indexfinger_apple > 0
            precision_grip_error = [precision_grip_error,i];
        elseif thumb_apple(3) > 0
            precision_grip_error = [precision_grip_error,i];
        end
    else
        precision_grip_error = [precision_grip_error,i];
    end
end

precision_grip_error = trial_start(unique(precision_grip_error));
error_rate3 = length(precision_grip_error)/trial_count;

end