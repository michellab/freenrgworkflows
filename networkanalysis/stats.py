#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Antonia Mey", "Julien Michel"
__email__ = "antonia.mey@ed.ac.uk"

import numpy as np
import networkx as nx
import scipy.stats
import copy

class freeEnergyStats(object):
    """docstring for freeEnergyStats"""
    def __init__(self, arg):
        super(freeEnergyStats, self).__init__()
        self.arg = arg
        

    def calculate_pi(self, series1, series2):
        sumwijcij = 0.0
        sumwij = 0.0

        keys = series1.keys()
        keys.sort()

        for i in range(0,len(keys)):
            keyi = keys[i]
            for j in range(i+1,len(keys)):
                keyj = keys[j]
                wij = abs(series1[keyj][0] - series1[keyi][0] )
                # print "series0 j %s series 0 i %s wij %s i %s j %s" % (series[0][j],series[0][i],wij,i,j)
                num =  (series1[keyj][0] - series1[keyi][0])
                den =  (series2[keyj][0] - series2[keyi][0] )
                #if den < 0.0001:
                #    den = 0.001
                #print num, den
                val = num / den
                # print val,serie[j],serie[i]
                if val > 0:
                    cij = 1.0
                elif val < 0:
                    cij = -1.0
                # print cij
                sumwijcij += wij*cij
                sumwij += wij
                # print i,j,series[0][j],serie[j],series[0][i],serie[i],val,wij*cij,wij
                # sys.exit(-1)
        PI = sumwijcij/sumwij
        #print PI
        return PI

    def calculate_r2 (self, series1, series2):
        r_value,p = scipy.stats.pearsonr(series1,series2)

        return r_value**2, r_value

    def calculate_tau(self, series1, series2):
        tau = scipy.stats.kendalltau(series1, series2)
        return tau[0]

    def calculate_mue(self, series1, series2 ):

        sumdev = 0.0
        for x in range(0,len(series1)):
            sumdev += abs( series1[x] - series2[x] )
        sumdev /= len(series1)

        #print sumdev
        return sumdev

    def perturb(self, data):
        r""" generates new set of data based on gauss distribution
        Parameters
        ----------
        data : nd.array(shape(datapoints,2))
            first column holding actual data, second error on data

        """
        repeat = np.zeros(np.shape(data))

        count = 0
        for d in data:
            val = d[0]
            err = d[1]
            if err != 0.0:
                val2 = np.random.normal(val, err)
            else:
                val2 = val
            repeat[count][0] = val2
            repeat[count][1] = err
            count = count + 1

        return repeat