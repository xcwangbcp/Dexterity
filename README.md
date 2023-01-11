# Dexterity
本项目是根据已经追踪好的2维（x，y）视频，通过手和苹果的两位坐标，判读被试抓取苹果用的时间，和错误率，苹果掉落率，其中错误里参考日本Nature的工作，是三种错误之和，slithit,wander,
gsp，由于我们装置的限制gsp因为不准确，没有更好的方法得到这个参数，后来放弃。同时由于现实的原因，苹果到边缘的距离也并不十分consistent，所以

Dexterity is the main funtion, which need to pepare the data into a folder, but data should be in the same month, 对于11号猴子，it will return e.g '11-08-09.mat',
’11-08-30.mat‘，每个mat文件了储存了那次实验的分析结果，同一只猴子的当月数据被平均了，具体数据存在了.xls里。
changeFilesName.m 可以批量修改文件的后缀，或者文件的名字等。
MutifileRead,可以批处理所有数据，e.g. 把2021年和2022年数据放到一个文件夹，输出的是同一个group之下的所有猴子在failIndex，speedIndex和dropIndex的平均值，并同时save到一个.xls
文件和.文件中。

统计部分：使用软件Jasp，如果操作可以参考我的石墨文档，需要一个csv文件。
