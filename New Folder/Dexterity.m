close all
clear 

filename       = uigetfile('*.csv','Pick a csv file to load in ');
[raw,Txt,RAW]  = xlsread(filename);
%1 base_x	base_y	base_p	
%2 MCP1_x	MCP1_y	MCP1_p	
%3 PIP1_x	PIP1_y	PIP1_p
%4 tip1_x	tip1_y	tip1_p
%5 MCP2_x	MCP2_y	MCP2_p
%6 PIP2_x	PIP2_y	PIP2_p
%7 DIP2_x	DIP2_y	DIP2_p
%8 tip2_x	tip2_y	tip2_p
filename_raw_hand = 'D:\Code\Dexterity\137-hand.csv';
filename_raw_apple = 'D:\Code\Dexterity\137-candy.csv';
Acrylic_Edge        = 'D:\Code\Dexterity\137LS_trimmed_cam1.csv';
raw_hand  = table2array(readtable(filename_raw_hand));
raw_apple = table2array(readtable(filename_raw_apple));
Edge_raw       = table2array(readtable(Acrylic_Edge));
Edge(:,1)  =  mean(Edge_raw(:,2));
Edge(:,2)  =  mean(Edge_raw(:,5));
Edge       =  mean(Edge);

M.base   = findcolum(Txt,'base_x','base_y','base_p');
M.MCP1   = findcolum(Txt,'MCP1_x','MCP1_y','MCP1_p');
M.PIP1   = findcolum(Txt,'PIP1_x','PIP1_y','PIP1_p');
M.tip1   = findcolum(Txt,'tip1_x','tip1_y','tip1_p');
M.MCP2   = findcolum(Txt,'MCP2_x','MCP2_y','MCP2_p');
M.PIP2   = findcolum(Txt,'PIP2_x','PIP2_y','PIP2_p');
M.DIP2   = findcolum(Txt,'DIP2_x','DIP2_y','DIP2_p');
M.tip2   = findcolum(Txt,'tip2_x','tip2_y','tip2_p');
index_tip= raw(:,M.tip2(1));
index_tip= movmean(index_tip,9);
touch = 0;time_fwd=[]; time_back=[];
for i=1:length(index_tip)-1
    if index_tip(i,1)<=Edge&&Edge<=index_tip(i+1,1) 
        
        touch = touch+1;
        time_fwd =[time_fwd,i];
    end
    if  index_tip(i,1)>=Edge&&Edge>=index_tip(i+1,1)
        time_back=[time_back,i];
    end
    
end

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