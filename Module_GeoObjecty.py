#!/usr/bin/env python
#-*- coding=utf-8 -*-
import time as Timetool
class GeoPoint:
    #基本的属性
    id=""
    x = 0
    y = 0
    isPickup = None
    isDropoff = None
    status = None
    time = "" #时间戳

    def setID(self,id):
        self.id=id
    def setX(self,x):
        self.x=float(x)-0.0025
    def setY(self,y):
        self.y=float(y)+0.0027
    def setTime(self,time):
#        self.time= Timetool.mktime(Timetool.strptime(time,'%Y/%m/%d %H:%M:%S'))
         self.time = time
    def setStatus(self,status):
        self.status=status
    def setIspickup(self,ispickup):
        self.isPickup = ispickup

    def __init__(self,point):
        self.setID(point[0])
        self.setY(point[1])
        self.setX(point[2])
        self.setStatus(point[3])
        self.setTime(point[4])

class Segment:
    id = ""
    start_time=""
    end_time=""
    length=""
    pickupnumber=0

