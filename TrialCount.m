function [trial_count] = TrialCount(raw_hand,raw_apple,raw_slit)
%hand from left only
%all plot() only for giving explanations
fps = 60;
raw_hand(find(raw_hand(:,24)<0.95),22:24) = nan;
raw_apple(find(raw_apple(:,3)<0.95),1:3) = nan;

raw_slit(find(raw_slit(:,4)<0.95),2:4) = nan;
raw_slit(find(raw_slit(:,7)<0.95),5:7) = nan;

raw_slit(find(raw_slit(:,10)<0.95),8:10) = nan;
raw_slit(find(raw_slit(:,13)<0.95),11:13) = nan;

left_top = [nanmean(raw_slit(:,2)),nanmean(raw_slit(:,3))];
left_bottom = [nanmean(raw_slit(:,5)),nanmean(raw_slit(:,6))];

TIP2 = movmean(raw_hand(:,22:23),9);
apple_start = raw_apple(:,1:2);
time = 1:size(apple_start,1);

%cleaning test
%nanmean(apple_start(:,1))
%nanstd(apple_start(:,1))

%original apple position
plot(apple_start(:,1),time,"y.")
hold on
%data cleaning by x value
apple_std_x = movstd(apple_start(:,1),9);
for i = 1:length(apple_std_x)
    if apple_std_x(i)>1 || isnan(apple_std_x(i))
        apple_start(i,:) = nan;
    end
end
start = isoutlier(apple_start(:,1));
for i = 1:size(apple_start,1)
    if start(i) == 1
        apple_start(i,:) = nan;
    end
end
plot(apple_start(:,1),time,"b.")
hold on
%data cleaning by y value
apple_std_y = movstd(apple_start(:,2),9);
for i = 1:length(apple_std_y)
    if apple_std_y(i)>1 || isnan(apple_std_y(i))
        apple_start(i,:) = nan;
    end
end
start = isoutlier(apple_start(:,2));
for i = 1:size(apple_start,1)
    if start(i) == 1
        apple_start(i,:) = nan;
    end
end
plot(apple_start(:,1),time,"k.")
%the cleaning cause an 8 frame data lost of apple x position(instead by
%NaN) on the start and the end of an independent moving, though it is
%shorter than 0.2s
plot(TIP2(:,1),time)

%cleaning test
%nanmean(apple_start(:,1))
%nanstd(apple_start(:,1))

slit_line_x = 0.5*(left_top(1)+left_bottom(1));
plot([slit_line_x,slit_line_x],[0,length(time)])

trial_count = 0;

for i = 1:size(apple_start,1)-7
    if TIP2(i,1) < slit_line_x && TIP2(i+1,1) >= slit_line_x %x pass  
        if sum(isnan(apple_start(i:i+7,1))) < 8 %apple on position
            if TIP2(i,2) > left_top(2) && TIP2(i,2) < left_bottom(2) %y in range
                if isempty(find(TIP2(i+1:end,1) < slit_line_x,1)+i) == 0 && find(TIP2(i+1:end,1) < slit_line_x,1) > 0.25*fps %record an end
                    if sum(isnan(TIP2(i:find(TIP2(i+1:end,1) < slit_line_x,1)+i))) == 0 %no lost data 
                        trial_count = trial_count+1;
                    end
                end
            end
        end
    end
end
          
end