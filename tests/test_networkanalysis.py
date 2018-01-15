import pytest
from networkanalysis import * 

def test_populate_pert_graph():
    pG = PerturbationGraph()
    pG.populate_pert_graph('tests/io/graph.csv')
    assert (pG.graph != None)
