# Dexterity
本项目是根据已经追踪好的2维（x，y）视频，通过手和苹果的两维坐标，判读被试抓取苹果用的时间，和错误率，苹果掉落率，其中错误率参考日本Nature的工作，是三种错误之和，slithit,wander,
gsp，由于我们装置的限制gsp因为不准确，没有更好的方法得到这个参数，后来放弃。同时由于现实的原因，苹果到边缘的距离也并不十分consistent，所以使用原始数据除以了distance，做了归一化。

1. Dexterity is the main funtion, which need to pepare the data into a folder, but data should be in the same month, 返回数据举例：如果11号猴子在8月有两次实验，会得到类似的结果 '11-08-09.mat',’11-08-30.mat‘，每个mat文件了储存了那次实验的分析结果，并且同一只猴子的当月数据以table格式存在了.xls里。

2. changeFilesName.m 可以批量修改文件的后缀，或者文件的名字等。

3. MutifileRead,可以批处理很多数据，处理对象是#1中生成的xls文件， 可以把2021年和2022年数据放到一个文件夹，输出的是同一个group之下的所有猴子在failIndex，speedIndex和dropIndex的平均值，并同时save到一个.xls文件和.mat文件中。

4. 统计部分：使用软件Jasp，如果操作可以参考我的石墨文档，需要主备一个csv文件。
