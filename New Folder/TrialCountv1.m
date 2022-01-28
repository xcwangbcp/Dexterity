function [trial_count,t] = TrialCountv1(raw_hand,raw_apple,slit_line,Edge)
%hand from left only
%all plot() only for giving explanations
%slit_line = [a,b] x=ay+b

raw_hand(find(raw_hand(:,24)<0.95),22:24) = nan;
raw_apple(find(raw_apple(:,9)<0.95),:) = nan;

TIP2 = movmean(raw_hand(:,22:23),9);
apple_start = raw_apple(:,7:8);
time = 1:size(apple_start,1);

%cleaning test
%nanmean(apple_start(:,1))
%nanstd(apple_start(:,1))

%original apple position
plot(apple_start(:,1),time,"r.")
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
plot(apple_start(:,1),time,"g.")
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
plot(apple_start(:,1),time,"b.")
%the cleaning cause an 8 frame data lost of apple x position(instead by
%NaN) on the start and the end of an independent moving, though it is
%shorter than 0.2s
plot(TIP2(:,1),time)

stem(Edge,time(end))
%cleaning test
%nanmean(apple_start(:,1))
%nanstd(apple_start(:,1))

slit_line_x = [];
for i = 1:size(TIP2,1)
    slit_line_x(i) = slit_line(1)*TIP2(i,2)+slit_line(2);
end
trial_count = 0;
time_fwd=[]; time_back=[];




p=0;q=0;
for i = 9:size(apple_start,1)
    if TIP2(i,1) <= Edge && TIP2(i+1,1) >= Edge%&&TIP2(i+10,1)>Edge%&& sum(isnan(apple_start(i-8:i))) > 0
        time_fwd =[time_fwd,i];
        trial_count = trial_count+1;
        hold on
        p=p+1;
        text(Edge-15,i,[num2str(p),'\rightarrow'],'Color','red','FontSize',15)
        j=i;
        while 1
            
            j=j+1
            if j>size(apple_start,1)
                break
            end
            if TIP2(j,1) >= Edge && TIP2(j+1,1) <= Edge
                
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

end