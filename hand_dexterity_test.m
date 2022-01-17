close all
clear 
Fs = 60;
% filename_likely = 'D:\Code\HandDexterity\Data\137-hand.csv';
filename_raw_hand = 'D:\Code\Dexterity\137-hand.csv';
filename_raw_apple = 'D:\Code\Dexterity\137-candy.csv';
raw_hand = table2array(readtable(filename_raw_hand));
raw_apple = table2array(readtable(filename_raw_apple));

%1. trial_count
slit_line = [0,320];
[trial_count] = TrialCount(raw_hand,raw_apple,slit_line)
%%
%2. success rate
[success] = Success(trial_count,apple_movement,apple_end)
[success_rate] = SuccessRate(trial_count,success)
%2. movement time
[movement_time] = MovementTime(TIP1_time,TIP2_time,PIP1_time,PIP2_time,MCP1_time,MCP2_time,success)