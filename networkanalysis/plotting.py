#!/usr/bin/env python
# -*- coding: utf-8 -*-

#!/usr/bin/env python

# This file is part of freenrgworkflows.
#
# Copyright 2016,2017 Julien Michel Lab, University of Edinburgh (UK)
#
# freenrgworkflows is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = "Antonia Mey"
__email__ = "antonia.mey@ed.ac.uk"

import matplotlib.pylab as plt
import matplotlib
import seaborn as sns
import copy
sns.set_style("ticks")
sns.set_context("notebook", font_scale = 2)

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
        print (r1_weight)
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