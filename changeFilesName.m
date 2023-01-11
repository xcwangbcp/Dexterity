function changeFilesName(fileOldMode,fileNewMode)
close all
clc;
% this function can change the file type , also can add some houzhui to the
% files
% eg fileMode:map4/csv,HorA:hand/apple changeFilesName('map4','hand')
% to change a lot of files's name fastly
fileM = ['*.' fileOldMode];
files = dir(fileM);    % 这个是文件存放的绝对路径
len   = length(files);  % 获取当前文件的长度，这里是1000

for i = 1: len        % 开始循环
    oldname = files(i).name;                     % 获取当前mat文件的名字
    temp    = ['-' fileNewMode '.' fileOldMode]; % 将数字转换为字符串
    newname = strcat(oldname(1:end-5),temp);     % 利用strcat函数进行字符串连接
    command = ['rename' 32 oldname 32 newname];  % 使用命令进行重命名
    status  = dos(command);                      % 调用dos命令
    if status == 0
        disp([oldname, ' rename to ', newname]);
    else
        disp([oldname, 'rename failed'])
    end
end
end
