#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

#输入为标准输入stdin
for line in sys.stdin:
# f = open("2")
# for line in f:
    #一行就是一个轨迹点
    line = line.strip()
    #以默认空格分隔单词到words列表
    words = line.split(",")
    content=words[0] # 先暂存id，与后面再一起组成整条记录
    pf = 0
    for word in words:
        #输出轨迹，2,30.658170,104.078292,0,2014/8/3 13:14:05
        if (pf==0):
            key = word
            pf = pf + 1
        else:
            content = content + "," + word
    print '%s\t%s' % (key,content)
# f.close()