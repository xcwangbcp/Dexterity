close all
clear 
Fs = 60;
% filename_likely = 'D:\Code\HandDexterity\Data\137-hand.csv';
filename_raw_hand = 'D:\Code\Dexterity\137-hand.csv';
filename_raw_apple = 'D:\Code\Dexterity\137-candy.csv';
Acrylic_Edge        = 'D:\Code\Dexterity\137LS_trimmed_cam1.csv';
raw_hand  = table2array(readtable(filename_raw_hand));
raw_apple = table2array(readtable(filename_raw_apple));
Edge_raw       = table2array(readtable(Acrylic_Edge));
Edge(:,1)  = mean(Edge_raw(:,2));
Edge(:,2)  = mean(Edge_raw(:,5));
Edge       = mean(Edge);

%1. trial_count
slit_line = [0,328];
[trial_count] = TrialCount(raw_hand,raw_apple,slit_line,Edge)
%%
%2. success rate
[success] = Success(trial_count,apple_movement,apple_end)
[success_rate] = SuccessRate(trial_count,success)
%2. movement time
[movement_time] = MovementTime(TIP1_time,TIP2_time,PIP1_time,PIP2_time,MCP1_time,MCP2_time,success)