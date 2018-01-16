import pytest
import warnings
import networkx as nx
from networkanalysis.networkanalysis import *

@pytest.fixture
def pG():
    return PerturbationGraph()

def test_populate_pert_graph_not_None(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert (pG.graph != None)

def test_compountList(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert ('FXR17' in pG.compoundList)

def test_double_call_pert_graph_not_None(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert(pG.populate_pert_graph('tests/io/graph.csv')==1)

#def test_symmetrize_graph():
#    newGraph = nx.read_edgelist('tests/io/graph.csv', delimiter=',', comments='#', nodetype=str, data=(('weight', float),('error',float)))

#def test_cycle_closure(pG, capsys):


#def test_compute_weighted_avg_paths():

#def test_compute_average_paths():

#def test_remove_compound_from_graph():

#def test_format_free_energies():
