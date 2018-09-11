#!/usr/bin/env python
# -*-coding=utf-8-*-
"""
-------------------------------------------------------------------------------
# Name:        mapMatcher
# Purpose:      This python script allows map matching (matching of track points to a network)
#               in arcpy using a Hidden Markov model with
#               probabilities parameterized based on spatial + network distances.
#               Follows the ideas in Newson, Krumm (2009):
#               "Hidden markov Map Matching through noise and sparseness"
#
#               Example usage under '__main__'
#
# Author:      Zhangjian follows Simon Scheider
#
# Created:     09/03/2018
# Copyright:   (c) zhangjian
# Licence:     <your licence>

The code is written in Python 2.7 and depends on:

* arcpy (ships with ArcGIS and its own Python 2.7)
* networkx (# python pip install networkx (https://networkx.github.io))
    (note: requires installing GDAL first, which can be obtained as a wheel from
    http://www.lfd.uci.edu/~gohlke/pythonlibs/ and then installed with pip locally:
    python pip install GDAL-2.1.3-cp27-cp27m-win32.whl
    )

#-------------------------------------------------------------------------------
"""

__author__ = "Zhangjian follows Simon Scheider"
__copyright__ = ""

import sys

try:
    from math import exp, sqrt
    import os
    from Module_GeoObjecty import Segment

#    import arcpy
    import shapefile
    import math
    from rtree import index
    import networkx as nx
    import time

except ImportError:
    #print "Error: missing one of the libraries (arcpy, networkx)"
    print("hello")
    sys.exit()

class MapMatchNArcpy:

    def __init__(self):
        self.idx = index.Index()
    def pointToSegment(self,x,y,x1,y1,x2,y2):
        cross = (x2 - x1) * (x - x1) + (y2 - y1) * (y - y1);
        if (cross <= 0):
            return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
        d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
        if (cross >= d2):
            return sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));
        r = cross / d2;
        px = x1 + (x2 - x1) * r;
        py = y1 + (y2 - y1) * r;
        return sqrt((x - px) * (x - px) + (py - y) * (py - y));

    def pointToSegments(self,point,segments):
        mindistance=float('inf')
        ii = 0
        for i in range(1,len(segments)):
            distance = self.pointToSegment(point.x,point.y,segments[ii][0],segments[ii][1],segments[i][0],segments[i][1])
            ii = i
            mindistance = distance if mindistance>distance else mindistance
        return mindistance

    def cleanPath(self, opt, endpoints):
        # removes redundant segments and segments that are unnecessary to form a path (crossings) in an iterative manner
        last = ()
        lastlast = ()
        optout = []
        for s in opt:
            if s != last:
                match = False
                if last != () and lastlast != ():
                    lastep = endpoints[last]
                    lastlastep = endpoints[lastlast]
                    sep = endpoints[s]
                    for j in lastlastep:
                        if lastep[0] == j:
                            for k in sep:
                                if lastep[1] == k:
                                    match = True
                        elif lastep[1] == j:
                            for k in sep:
                                if lastep[0] == k:
                                    match = True
                elif last != ():
                    sep = endpoints[s]
                    lastep = endpoints[last]
                    for k in sep:
                        if lastep[1] == k or lastep[0] == k:
                            match = True
                if match:
                    optout.append(last)
                if s == opt[-1]:
                    # print "add final segment:"+str(s)
                    optout.append(s)
            lastlast = last
            last = s
        # print "final length: "+str(len(optout))
        return optout

    def mapMatch(self,track, segments, decayconstantNet=30, decayConstantEu=10, maxDist=50, addfullpath=True):

        """
        The main method. Based on the Viterbi algorithm for Hidden Markov models,
        The main method. Based on the Viterbi algorithm for Hidden Markov models,
        see https://en.wikipedia.org/wiki/Viterbi_algorithm.
        It gets trackpoints and segments, and returns the most probable segment path (a list of segments) for the list of points.
        Inputs:
            @param track = a shape file (filename) representing a track, can also be unprojected (WGS84)
            @param segments = a shape file of network segments, should be projected (in meter) to compute Euclidean distances properly (e.g. GCS Amersfoord)
            @param decayconstantNet (optional) = the network distance (in meter) after which the match probability falls under 0.34 (exponential decay). (note this is the inverse of lambda).
            This depends on the point frequency of the track (how far are track points separated?)
            @param decayConstantEu (optional) = the Euclidean distance (in meter) after which the match probability falls under 0.34 (exponential decay). (note this is the inverse of lambda).
            This depends on the positional error of the track points (how far can points deviate from their true position?)
            @param maxDist (optional) = the Euclidean distance threshold (in meter) for taking into account segments candidates.
            @param addfullpath (optional, True or False) = whether a contiguous full segment path should be outputted. If not, a 1-to-1 list of segments matching each track point is outputted.

        note: depending on the type of movement, optional parameters need to be fine tuned to get optimal results.
        """
        # Make sure passed in parameters are floats
        decayconstantNet = float(decayconstantNet)
        decayConstantEu = float(decayConstantEu)
        maxDist = float(maxDist)

        # gest start time
        start_time = time.time()

        # this array stores, for each point in a track, probability distributions over segments, together with the (most probable) predecessor segment taking into account a network distance
        V = [{}]

        # get track points, build network graph (graph, endpoints, lengths) and get segment info from arcpy
        points = self.getTrackPoints(track)
        r = self.getSegmentInfo(segments)
        endpoints = r[0]
        lengths = r[1]
        graph = self.getNetworkGraph(segments, lengths)
        pathnodes = []  # set of pathnodes to prevent loops

        # init first point
        sc = self.getSegmentCandidates(points[0], segments, decayConstantEu, maxDist)
        for s in sc:
            V[0][s] = {"prob": sc[s], "prev": None, "path": [], "pathnodes": []}
        # Run Viterbi when t > 0
        for t in range(1, len(points)):
            V.append({})
            # Store previous segment candidates
            lastsc = sc
            # Get segment candidates and their a-priori probabilities (based on Euclidean distance for current point t)
            sc = self.getSegmentCandidates(points[t], segments, decayConstantEu, maxDist)
            if len(sc)<=0:
                return None
            for s in sc:
                max_tr_prob = 0
                prev_ss = None
                path = []
                for prev_s in lastsc:
                    # determine the highest network transition probability from previous candidates to s and get the corresponding network path
                    pathnodes = V[t - 1][prev_s]["pathnodes"][-10:]
                    n = self.getNetworkTransP(prev_s, s, graph, endpoints, lengths, pathnodes, decayconstantNet)
                    np = n[0]  # This is the network transition probability
                    tr_prob = V[t - 1][prev_s]["prob"] * np
                    # this selects the most probable predecessor candidate and the path to it
                    if tr_prob > max_tr_prob:
                        max_tr_prob = tr_prob
                        prev_ss = prev_s
                        path = n[1]
                        if n[2] != None:
                            pathnodes.append(n[2])
                # The final probability of a candidate is the product of a-priori and network transitional probability
                max_prob = sc[s] * max_tr_prob
                V[t][s] = {"prob": max_prob, "prev": prev_ss, "path": path, "pathnodes": pathnodes}

            # Now max standardize all p-values to prevent running out of digits
            maxv = max(value["prob"] for value in V[t].values())
            maxv = (1 if maxv == 0 else maxv)
            for s in V[t].keys():
                V[t][s]["prob"] = V[t][s]["prob"] / maxv

        intertime1 = time.time()
        print("--- Viterbi forward: %s seconds ---" % (intertime1 - start_time))
        # print V

        # opt is the result: a list of (matched) segments [s1, s2, s3,...] in the exact order of the point track: [p1, p2, p3,...]
        opt = []

        # get the highest probability at the end of the track
        if len(V[-1].values())==0 or V[-1]==None:
            return None
        max_prob = max(value["prob"] for value in V[-1].values())
        previous = None
        if max_prob == 0:
            print " probabilities fall to zero (network distances in data are too large, try increasing network decay parameter)"
            return None  # if prob=0,then end this regulation.
        # Get most probable ending state and its backtrack
        # recode the maprelation of point and route segment
        PR = []
        for st, data in V[-1].items():
            if data["prob"] == max_prob:
                opt.append(st)
                previous = st
                PR.append(str(st))
                break
        ##    print  " previous: "+str(previous)
        ##    print  " max_prob: "+str(max_prob)
        ##    print  " V -1: "+str(V[-1].items())

        # Follow the backtrack till the first observation to fish out most probable states and corresponding paths

        for t in range(len(V) - 2, -1, -1):
            # Get the subpath between last and most probable previous segment and add it to the resulting path
            if V[t + 1].has_key(previous):  # if can match then go ,else do nothing
                path = V[t + 1][previous]["path"]
                opt[0:0] = (path if path != None else [])
                # Insert the previous segment
                opt.insert(0, V[t + 1][previous]["prev"])
                PR.insert(0, str(previous))
                previous = V[t + 1][previous]["prev"]

        intertime2 = time.time()
        print("--- Viterbi backtracking: %s seconds ---" % (intertime2 - intertime1))

        # Clean the path (remove double segments and crossings) (only in full path option)
        print "path length before cleaning :" + str(len(opt))
        opt = self.cleanPath(opt, endpoints)
        intertime3 = time.time()
        print("--- Path cleaning: %s seconds ---" % (intertime3 - intertime2))
        print "final length: " + str(len(opt))
        pointstr = [str(g.x) + ' ' + str(g.y) for g in points]
        optstr = [str(i) for i in opt]
        print 'The path for points [' + ' '.join(pointstr) + '] is: '
        print '[' + ' '.join(optstr) + '] with highest probability of %s' % max_prob
        return PR

    def getPDProbability(self,dist, decayconstant=10):
        """
        The probability that given a certain distance between points and segments, the point is on the segment
        This needs to be parameterized
        Turn difference into a probability with exponential decay function
        """
        decayconstant = float(decayconstant)
        dist = float(dist)
        try:
            p = 1 if dist == 0 else round(1 / exp(dist / decayconstant), 4)
        except OverflowError:
            p = round(1 / float('inf'), 2)
        return p


    def getSegmentCandidates(self,point, segments, decayConstantEu, maxdist=50):
        """
        Returns closest segment candidates with a-priori probabilities.
        Based on maximal spatial distance of segments from point.
        """
        # print "Neighbors of point "+str(p.X) +' '+ str(p.Y)+" : "
        # Select all segments within max distance
        candidates = {}
        left = point.x-maxdist
        bottom = point.y-maxdist
        right = point.x+ maxdist
        top =  point.y+ maxdist
        candidateSegmentObjects = [n for n in self.idx.intersection((left, bottom, right, top), objects=True)]
        for SegmentObjects in candidateSegmentObjects:
            dist = self.pointToSegments(point,SegmentObjects.object)
            candidates[SegmentObjects.id] = self.getPDProbability(dist, decayConstantEu)
        return candidates

    def getNDProbability(self,dist, decayconstant=30):
        """
        The probability that given a certain network distance between segments, one is the successor of the other in a track
        This needs to be parameterized
        Turn difference into a probability  with exponential decay function
        """
        decayconstant = float(decayconstant)
        dist = float(dist)
        try:
            p = 1 if dist == 0 else round(1 / exp(dist / decayconstant), 2)
        except OverflowError:
            p = round(1 / float('inf'), 2)
        return p


    def getNetworkTransP(self,s1, s2, graph, endpoints, segmentlengths, pathnodes, decayconstantNet):
        """
        Returns transition probability of going from segment s1 to s2, based on network distance of segments, as well as corresponding path
        """
        subpath = []
        s1_point = None
        s2_point = None

        if s1 == s2:
            dist = 0
        else:
            # Obtain edges (tuples of endpoints) for segment identifiers
            s1_edge = endpoints[s1]
            s2_edge = endpoints[s2]

            s1_l = segmentlengths[s1]
            s2_l = segmentlengths[s2]

            # This determines segment endpoints of the two segments that are closest to each other
            minpair = [0, 0, 100000]
            for i in range(0, 2):
                for j in range(0, 2):
                    d = round(self.pointdistance(s1_edge[i], s2_edge[j]), 5)
                    if d < minpair[2]:
                        minpair = [i, j, d]
            s1_point = s1_edge[minpair[0]]
            s2_point = s2_edge[minpair[1]]

            ##        if (s1_point in pathnodes or s2_point in pathnodes): # Avoid paths reusing an old pathnode (to prevent self-crossings)
            ##            dist = 100
            ##        else:
            if s1_point == s2_point:
                # If segments are touching, use a small network distance
                dist = 0.00005
            else:
                try:
                    # Compute a shortest path (using segment length) on graph where segment endpoints are nodes and segments are (undirected) edges
                    if graph.has_node(s1_point) and graph.has_node(s2_point):
                        dist = nx.shortest_path_length(graph, s1_point, s2_point, weight='length')
                        path = nx.shortest_path(graph, s1_point, s2_point, weight='length')
                        # get path edges
                        path_edges = zip(path, path[1:])
                        # print "edges: "+str(path_edges)
                        subpath = []
                        # get object ids for path edges
                        for e in path_edges:
                            oid = graph[e[0]][e[1]]["OBJECTID"]
                            subpath.append(oid)
                        # print "oid path:"+str(subpath)
                    else:
                        # print "node not in segment graph!"
                        dist = 3 * decayconstantNet  # 600
                except nx.NetworkXNoPath:
                    # print 'no path available, assume a large distance'
                    dist = 3 * decayconstantNet  # 700
        # print "network distance between "+str(s1) + ' and '+ str(s2) + ' = '+str(dist)
        return (self.getNDProbability(dist, decayconstantNet), subpath, s2_point)


    def pointdistance(self,p1, p2):
        # This Eucl distance can only be used for projected coordinate systems
        dist = sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
        return dist


    def getTrackPoints(self,track):
        """
        Turns track shapefile into a list of point geometries, reprojecting to the planar RS of the network file
        """
        trackpoints = []
        for point in track:
            trackpoints.append(point)
        return trackpoints


    def getNetworkGraph(self,segmentsfile, segmentlengths):
        """
        Builds a networkx graph from the network file, inluding segment length taken from arcpy.
        It selects the largest connected component of the network (to prevent errors from routing between unconnected parts)
        """
        # generate the full network path for GDAL to be able to read the file
        if (os.path.exists(segmentsfile)):
            g = nx.read_shp(segmentsfile)
            # This selects the largest connected component of the graph
            sg = list(nx.connected_component_subgraphs(g.to_undirected()))[0]
            print "graph size (excluding unconnected parts): " + str(len(g))
            # Get the length for each road segment and append it as an attribute to the edges in the graph.
            for n0, n1 in sg.edges():
                oid = sg[n0][n1]["OBJECTID"]
                sg[n0][n1]['length'] = segmentlengths[oid]
            return sg
        else:
            print "network file not found on path: " + segmentsfile

    def lengOfSegment(self,Points):
        ii=0
        sumlenth = 0
        for i in range(1,len(Points)):
            sumlenth = sumlenth + self.pointdistance(Points[ii],Points[i])
            ii = i
        return sumlenth


    def getSegmentInfo(self,segmentsfile):
        """
        Builds a dictionary for looking up endpoints of network segments (needed only because networkx graph identifies edges by nodes)
        """
        # construct rtree instance
        if (os.path.exists(segmentsfile)):

            sf = shapefile.Reader(segmentsfile)
            shapeRecs = sf.shapeRecords()
            endpoints = {}
            segmentlengths = {}
            for segment in shapeRecs:
                id = segment.record[0]
                points = segment.shape.points
                bbox = segment.shape.bbox
                self.idx.insert(id,bbox,points)
                endpoints[id] = ((points[0][0],points[0][1]), (points[-1][0], points[-1][1]))
                segmentlengths[id] = self.lengOfSegment(points)
            print "Number of segments: " + str(len(endpoints))
            return (endpoints, segmentlengths)
        else:
            print "segment file does not exist!"

    if __name__ == '__main__':
       sf = shapefile.Reader("D:\\shp\\sichuan\\mapmatch\\chengdu.shp")
       shapeRecs = sf.shapeRecords()
       print (pointToSegment(2,1,0,0,1,0))
    #
