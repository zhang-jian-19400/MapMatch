#!/usr/bin/env python
# -*- coding:utf-8 -*-

from operator import itemgetter
import sys
from Module_GeoObjecty import GeoPoint
from Module_GeoObjecty import Segment
import MapMatchNoArcpy
#import arcpy
current_key = None
key = None
trajectory = []

def divideSegment(trajectory):
    segments = []
    segment = []
    alpha = 1200 #超过半个小时的时间间隔，则分段
    lastpoint = trajectory[0]
    for point in trajectory:
        if (int(point.time) -int(lastpoint.time)>alpha) or (str(point.status) != str(lastpoint.status)):
            if (str(point.status) != str(lastpoint.status) and (str(point.status) == "1")):
                point.setIspickup("1")
            else:
                point.setIspickup("0")
            segments.append(segment)
            segment=[]
            segment.append(point)
        else:
            point.setIspickup("0")
            segment.append(point)
        lastpoint = point
    return segments

def countSegment(points,segments):
    ResultS = []
    lastsegment = segments[0]
    m = MapMatchNoArcpy.MapMatchNArcpy()
    track=[]
    number = 0
    pointset = []
    for i in range(0,len(points)):
        currentseg = segments[i]
        if (currentseg==lastsegment):
           if points[i].isPickup == "1":
              number = number + 1
           track.append(points[i])
           pointset.append((points[i].x,points[i].y))
        elif(len(track)>0):
            # record some segment info
            s = Segment()
            s.id = lastsegment
            s.start_time = track[0].time
            s.end_time = track[-1].time
            s.length = m.lengOfSegment(pointset) #may cause 0 length.
            s.pickupnumber = number
            ResultS.append(s)
            # clean and record current point
            pointset=[]
            track = []
            number = 0
            track.append(points[i])
            pointset.append((points[i].x, points[i].y))
        lastsegment = currentseg
    s = Segment()
    s.start_time = track[0].time
    s.end_time = track[-1].time
    s.length = m.lengOfSegment(pointset)
    s.pickupnumber = number
    s.id = lastsegment
    ResultS.append(s)
    return ResultS

def prn_obj(objs):
    content = ""
    for obj in objs:
        content = (content + ';'.join(['%s:%s' % item for item in obj.__dict__.items()]) + '\n')
    return content

def mapmatch(segments):
#    arcpy.env.workspace = 'D:\\shp\\sichuan\\mapmatch'
    content = ""
    for points in segments:
        print("段的点的个数为 %s") % len(points)
        m = MapMatchNoArcpy.MapMatchNArcpy()
        #path = m.mapMatch(points,'D:\\shp\\sichuan\\mapmatch\\chengdu.shp', 0.0003, 0.0003, 0.001)
        #path =m.mapMatch(points,'D:\\shp\\sichuan\\mapmatch\\chengdu.shp', 0.0003, 0.0003, 0.001)
        path = m.mapMatch(points, '/home/hadoop/zj/MapMatchParallel/chengdu.shp', 0.0003, 0.0003, 0.001)
        if path==None:
            print("no match")
        else:
            content =content + prn_obj(countSegment(points,path))
    return content
#获取标准输入，即mapper.py的标准输出
print("开始读轨迹信息。。。")
number = 0
for line in sys.stdin:
    number = number + 1
    if number==500:
        print("读到了前500行。。。")
# f = open("all")
# flag=False
# for line in f:
    #删除开头和结尾的空行
    line = line.strip()
    #解析mapper.py输出作为程序的输入，以tab作为分隔符
    key,content = line.split(',',1)

    #转换content字符串到GeoPoint对象
    # if key=='10055':
    #     flag = True
    #     continue
    #
    # if flag !=True:
    #     continue
    try:
        content = str(content)
        point =  GeoPoint(content.split(",")[0:5])
    except Exception:
        continue

#要求mapper.py的输出做排序（sort）操作，以便对连续的word做判断
    if current_key == key:
        #就把该记录记好
        trajectory.append(point)
    else:
        if current_key:
            #就要把当前组中的元素处理完，并输出。
            # 并按照时间排好序
            trajectory = sorted(trajectory, key=lambda item: item.time)
            #其处理工作包括：分段。轨迹地图映射。并统计出相应轨迹段的{轨迹段：{1：{时间：通行时间；时间：通行时间}:2：{时间} } }
            print("开始分段。。。")
            segments = divideSegment(trajectory)
            print '%s\t%s' % (current_key, str(mapmatch(segments)))
            # info = {current_key : str(mapmatch(segments))}
            # f = open('./result.txt', 'a+')
            # f.write(str(info))
            # f.close()
        #重新开始记录相关信息
        del trajectory
        trajectory = []
        trajectory.append(point)
        current_key = key

#计算最后一个人的轨迹
if current_key == key:
    # 就要把当前组中的元素处理完，并输出。
    # 并按照时间排好序
    sorted(trajectory, key=lambda item: item.time)
    # 其处理工作包括：分段。轨迹地图映射。并统计出相应轨迹段的{轨迹段：{1：{时间：通行时间；时间：通行时间}:2：{时间} } }
    segments = divideSegment(trajectory)
    print '%s\t%s' % (current_key, str(mapmatch(segments)))
    # info = {current_key: str(mapmatch(segments))}
#     # f = open('./result.txt', 'a+')
#     # f.write(str(info))
#     # f.close()
print("处理完成！")
#f.close()
