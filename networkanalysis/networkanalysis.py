#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Antonia Mey"
__email__ = "antonia.mey@ed.ac.uk"


import numpy as np
import networkx as nx
import matplotlib.pylab as plt
import matplotlib
import seaborn as sns
import scipy.stats
import copy
sns.set_style("ticks")
sns.set_context("notebook", font_scale = 2)

class PertubrationGraph(object):
    """Populates a directed free energy perturbation graph"""
    def __init__(self, ):
        self._graph = None
        self._pathAverages = []
        self._weightedPathAverages = []


    def populate_pert_graph(self, filename, delimiter=',', comments='#', nodetype=str, data=(('weight',float),('error',float))):
        r"""
        Parameters
        ----------
        filename : String
            filename of the forward and backward perturbation generated from simulation output
            Filestructure should be:
            node1,node2,DG,eDG,other_attributes
        delimiter : String
            delimiter for network file 
            Default = ','
        comments : String
            Symbol used for comments in network file
            Default = '#'
        nodetype : String
            All nodes are usually identified by the compound name
        data : list
            Default, weight and error on Free energies of node
        """
        self._graph = nx.read_edgelist(filename, delimiter=delimiter, comments=comments, create_using=nx.DiGraph(),nodetype=nodetype, data=data)


    def symmetrize_graph(self, graph):
        r"""
        Parameters
        ----------
        graph : networkx graph
            directed networkx graph

        Returns
        -------
        graph : networkx graph
            returns directed graph where, if not both a forward and backward edge are present a symmetrized reverse edge is included
        """
        for u,v,w in graph.edges(data=True):
            if not graph.has_edge(v,u):
                assymetric_w = -w['weight']
                assymetric_e = -w['error']
                graph.add_edge(v,u,weight=assymetric_w, error = assymetric_e)
        return graph

    def save_average_paths(self, filename, pathAverages):
        f = open(filename, 'w')
        for d in pathAverages:
            for k,v in d.iteritems():
                if k == 'error':
                    error = v
                else:
                    r_energy_k = k
                    r_energy_v = v
            f.write(str(r_energy_k)+', '+str(r_energy_v)+', '+str(error)+'\n')
        f.close()

    def compute_average_paths(self, target_node):
        r"""
        Parameters
        ----------
        target_node : string
            node to which all possible paths are computed
        """
        #Get all relative free energies with respect to node x
        for n in self._graph.nodes():
            paths = nx.shortest_simple_paths(self._graph,target_node , n)
            err_list = []
            sum_list = []
            for p in paths:
                sum = 0
                error = 0.0
                for node in range(len(p)-1):
                    sum= sum+ self._graph.get_edge_data(p[node], p[node+1])['weight']
                    error = error+self._graph.get_edge_data(p[node], p[node+1])['error']**2
                sum_list.append(sum)
                err_list.append(error)
                error = np.sqrt(error)
                #print r'DDG for path %s is %f +/- %f kcal/mol' %(p, sum, error)
            avg_sum = np.mean(np.array(sum_list))
            avg_err = np.mean(np.array(err_list))
            avg_std = np.std(np.array(sum_list))
            #print ("Average sum for path to %s is %f " %(n,avg_sum))
            a = {str(n):avg_sum}
            a['error']=avg_std
            #a['error']=sqrt(avg_err)
            self._pathAverages.append(a)

    def compute_weighted_avg_paths(self, target_node):
        #Get all relative free energies with respect to node x
        a = {target_node:0.0}
        a['error']=0.0
        self._weightedPathAverages.append(a)
        for n in self._graph.nodes():
            if n == target_node:
                continue
            #print "=============================="
            #print "Path: "+str(n)+"~"+str(target_node)
            paths = list(nx.shortest_simple_paths(self._graph,target_node , n))
            #print "There are "+str(len(paths))+" simple paths."
            err_list = []
            sum_list = []
            for p in paths:
                summing = 0
                error = 0.0
                for node in range(len(p)-1):
                    summing= summing+ self._graph.get_edge_data(p[node], p[node+1])['weight']
                    error = error+self._graph.get_edge_data(p[node], p[node+1])['error']**2
                sum_list.append(summing)
                error = np.sqrt(error)
                err_list.append(error)
                #print 'path sum is: ' +str(summing)
                #print 'path error is: '+str(error)
            err_list = np.array(err_list)
            sum_weights = np.sum(1.0/err_list)
            #print 'sum weights: '+str(sum_weights)
            #print 1.0/err_list
            path_weights = (1.0/err_list)/sum_weights
            #print path_weights
            avg_sum = 0.0
            avg_err = 0.0
            for i in range(len(sum_list)):
                s = sum_list[i]
                avg_sum = avg_sum +(path_weights[i]*s)
                avg_err = avg_err+path_weights[i]*err_list[i]**2
            avg_err = np.sqrt(avg_err)
            a = {str(n):avg_sum}
            a['error']=avg_err
            self._weightedPathAverages.append(a)
            #print "path average is: "+str(np.mean(sum_list))
            #print "path weighted average is: "+str(avg_sum)
            #print "sum err is: " +str(avg_err)
            #print "================"


    def get_cycles(self):
        #cycle closure
        cyc = nx.simple_cycles(self._graph)
        for c in cyc:
            sum = 0
            error = 0
            if len(c)>2:
                sum = self._graph.get_edge_data(c[-1], c[0])['weight']
                error = (self._graph.get_edge_data(c[-1], c[0])['error'])**2
                for node in range(len(c)-1):
                    sum= sum+ self._graph.get_edge_data(c[node], c[node+1])['weight']
                    error = error +(self._graph.get_edge_data(c[node], c[node+1])['error'])**2
                error = np.sqrt(error)
            print 'DDG for cycle %s is %.2f ± %.2f kcal/mol' %(c,sum,error)

    @property
    def graph(self):
        return self._graph

    @property
    def pathAverages(self):
        r"""
        Return
        ------
        pathAverages : dictionary
            dictionary containing averaged free energies for each path, with paths weighted in the same way 
        """
        return self._pathAverages

    @property
    def weightedPathAverages(self):
        return self._weightedPathAverages

class FreeEnergyPlotter(object):
    """docstring for FreeEnergyPlotter"""
    def __init__(self, arg):
        super(FreeEnergyPlotter, self).__init__()
        self.arg = arg
        

    def plot_bar_plot(self, set_ddg, graph_est, keys):
        r1_weight = []
        r2_weight = []
        labels = []
        for e in keys:
            for i in graph_ddg:
                if i.has_key(e):
                    r1_weight.append(i[e])
        for e in keys:
            for i in graph_est:
                if i.has_key(e):
                    r2_weight.append(i[e])
        labels = keys

        N = len(r1_weight)*2

        ind = np.arange(0,N,2)  # the x locations for the groups
        width = 0.35*2       # the width of the bars
        fig, ax = plt.subplots(figsize=(8,6))

        rects1 = ax.bar(ind, r1_weight, width, color=sns.xkcd_rgb["pale red"])
        rects2 = ax.bar(ind + width, r2_weight, width, color=sns.xkcd_rgb["denim blue"])

        # add some text for labels, title and axes ticks
        ax.set_ylabel(r'$\Delta \Delta G$ in [kcal/mol]', fontsize=15)
        ax.set_xticks(ind + width)
        ax.set_xticklabels(labels, fontsize=15)

        ax.legend((rects1[0], rects2[0]), ('experimental', 'computational'), fontsize=15)
        sns.despine()

    def plot_bar_plot_no_dic(self, r1_weight, r2_weight, keys):
        labels = keys

        N = len(r1_weight)*2

        ind = np.arange(0,N)  # the x locations for the groups
        width = 0.35*2       # the width of the bars
        fig, ax = plt.subplots(figsize=(7,7))

        rects1 = ax.bar(ind, r1_weight[:,0], width, yerr = r1_weight[:,1], color=sns.xkcd_rgb["pale red"])
        rects2 = ax.bar(ind + width, r2_weight[:,0], width, yerr = r2_weight[:,1], color=sns.xkcd_rgb["denim blue"])

        # add some text for labels, title and axes ticks
        ax.set_ylabel(r'$\Delta \Delta G$ in [kcal/mol]', fontsize=20)
        ax.set_xticks(ind + width)
        ax.set_xticklabels(labels, fontsize=20)

        #ax.legend((rects1[0], rects2[0]), ('experimental', 'computational'), fontsize=20)
        sns.despine()

    def plot_bar_plot_no_dic_graph2(self, r1_weight, r2_weight, r3_weight, keys):
        labels = keys

        N = len(r1_weight)*2

        ind = np.arange(0,N,2)  # the x locations for the groups
        width = 0.25*2       # the width of the bars
        fig, ax = plt.subplots(figsize=(7,7))
        print r1_weight
        data1 = r2_weight[:,0]
        data2 = r3_weight[:,0]
        rects1 = ax.bar(ind, r1_weight, width, color=sns.xkcd_rgb["pale red"])
        rects2 = ax.bar(ind + width, data1, width,yerr = r2_weight[:,1], color=sns.xkcd_rgb["denim blue"], error_kw={'ecolor':'black',    # error-bars colour
                              'linewidth':2})
        rects3 = ax.bar(ind + 2*width, data2, width, yerr = r3_weight[:,1], color='#759D70', error_kw={'ecolor':'black',    # error-bars colour
                              'linewidth':2})

        # add some text for labels, title and axes ticks
        ax.set_ylabel(r'$\Delta \Delta G$ in [kcal/mol]', fontsize=20)
        ax.set_xticks(ind + width+width*0.5)
        ax.set_xticklabels(labels, fontsize=20)

        #ax.legend((rects1[0], rects2[0]), ('experimental', 'computational'), fontsize=20)
        sns.despine()

    def plot_hist(self, edges, hist, label, color, alpha):
        halfwidth = (edges[1]-edges[0])/2
        centers= edges[:-1]+halfwidth
        fig = plt.figure(figsize=(8,4))
        plt.plot(centers, hist, color=color)
        plt.fill_between(centers, hist, facecolor=color, alpha=alpha)
        plt.xlabel(label)
        plt.ylabel('P(%s)'%label)
        sns.despine()

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


