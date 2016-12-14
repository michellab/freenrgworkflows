#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Antonia Mey"
__email__ = "antonia.mey@ed.ac.uk"


import numpy as np
import networkx as nx
import scipy.stats
import copy


class ExperimentalData(object):
    """docstring for ExperimentalData"""
    def __init__(self):
        self._DG_in_kcal = {}
        self._DG_in_kJ = {}
        self._ic50s = []
        self._referenceCompound = None
        self._kT = 0.0019872041*300

    def from_IC50s(self, filename, reference = None):
        r"""
        filename : string
            file containing ic50 data, format - compound name, ic50 value, error
        """
        f = open(filename, 'r')
        for line in f.readlines():
            curr_ic50 = {}
            fields = line.split(',')
            curr_ic50[fields[0]] = float(fields[1].strip())
            self._ic50s.append(curr_ic50)
            #TODO: check next entry, is this an error? Something else? 

        if reference is not None:
            self._referenceCompound = reference
        #for k in self._ic50s.keys():
        #    r = float(self._ic50s.get(k))/float(self._ic50s.get(self._referenceCompound))
        #    DDG = self.kT*np.log(r)
        print (self._ic50s)
        #for i in self._ic50s:



    def from_free_energies(self, filename, kcal=True):
        r"""
        filename : string
            file containing ic50 data, format - compound name, free energy_relative_to , error
        """
        #if kcal:
            #self._DG_in_kcal = np.loadtxt(filename)
        #else:
            #self._DG_in_kJ = np.loadtxt(filename)

    @property
    def ic50s(self):
        return self._ic50s

    @property
    def freeEnergiesInKcal(self):
        return self._DG_in_kcal

    @property
    def freeEnergiesInKJmol(self):
        return self._DG_in_kJ

