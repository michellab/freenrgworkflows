import pytest
import warnings
import networkx as nx
from networkanalysis.networkanalysis import *


@pytest.figure
def nA():
    return NetworkAnalyser()

def test_compoundList(nA):
    nA.read_perturbations('tests/io/graph.csv')
    assert('FXR17' in nA.compoundList)

def test_dG_simple(nA):
    nA.read_perturbations('tests/io/simple.csv')
    x,y = nA.dG()
    validate_results(x,[-0.5,0.5])
    validate_results(y,[0.4,0.4],0.05)

def test_perfectcycle(nA):
    nA.read_perturbations("tests/io/perfectcycle.csv")
    x,y = nA.dG()
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [0.4, 0.4, 0.4], 0.05)

def test_inconsistentcycle():
    nA.read_perturbations("tests/io/inconsistentcycle.csv")
    x,y = nA.dG()
    validate_results(x, [0, 0, 0])
    validate_results(y, [0.4, 0.4, 0.4], 0.05)

def test_inconsistentcycle_weights():
    nA.read_perturbations("tests/io/inconsistentcycle_weights.csv")
    x,y = nA.dG()
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [0.44, 0.44, 0.38], 0.05)

def test_large_hysteresis():
    nA.read_perturbations("tests/io/large_hysteresis.csv")
    x, y = nA.dG()
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [0.66, 0.66, 0.66], 0.05)

def test_vlarge_hysteresis():
    nA.read_perturbations("tests/io/vlarge_hysteresis.csv")
    x, y = nA.dG()
    validate_results(x, [-0.5, 0.5, 0])
    validate_results(y, [1.60, 1.60, 1.60], 0.05)

def test_large_cycle():
    nA.read_perturbations("tests/io/large_cycle.csv")
    x, y = nA.dG()
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0])
    validate_results(y, [0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65], 0.05)

def test_large_cycle_poor_link():
    nA.read_perturbations("tests/io/freenrg_large_cycle_poor_link.csv")
    x, y = nA.dG()
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0], 0.01)
    # These are the error values when the "poor" link is deleted: we should get the same when it
    # is present but has a very low weight
    validate_results(y, [1.1772, 0.9524, 0.7728, 0.6652, 0.6612, 0.7720, 0.9576, 1.1777], 0.05)

def test_DG_noise_handling():
    nA.read_perturbations("tests/io/noise0.csv")
    x, y = nA.dG()
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0, 1, 2, 3, 4, 4])
    # Errors: given the network, m4 should have the lowest error and m6 the highest
    validate_results(y, [0.54, 0.53, 0.43, 0.35, 0.53, 0.73], 0.05)

    # Same network, random error with std dev 0.5 added to all deltag measurements
    nA.read_perturbations("tests/io/noise0.5.csv")
    x, y = nA.dG()
    x = [t - x[0] for t in x]  # Set the first mol to zero
    validate_results(x, [0.0, 0.935, 2.425, 3.365, 4.64, 4.57])
    # Errors: given the network, m4 should have the lowest error and m6 the highest
    validate_results(y, [0.845, 0.902, 0.591, 0.564, 0.779, 0.864], 0.05)

def validate_results(x, xcorrect, delta=0.01):
    print("Checking",x,"against",xcorrect,"delta",delta)
    assert len(x) == len(xcorrect)
    for i in range(len(x)):
        assert x[i] == pytest.approx(xcorrect[i], abs=delta)

@pytest.fixture
def pG():
    return PerturbationGraph()


def test_populate_pert_graph_not_None(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert (pG.graph != None)


def test_compoundList(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    assert ('FXR17' in pG.compoundList)


def test_double_call_pert_graph_not_None(pG):
    pG.populate_pert_graph('tests/io/graph.csv')
    warn_string = 'Warning...........Use the method add_data_to_graph, to add further data to an existing graph'
    with pytest.warns(UserWarning) as warnmessg:
        pG.populate_pert_graph('tests/io/graph.csv')
    assert len(warnmessg) == 1
    assert warnmessg[0].message.args[0] == warn_string

def test_perfect_graph(pG):
    pG.populate_pert_graph('tests/io/test_perfect_graph.csv')
    G = pG.graph
    for edge in G.edges:
        data = G.get_edge_data(edge[0], edge[1])
        assert data['error'] is not 0.0
# def test_symmetrize_graph():
#    newGraph = nx.read_edgelist('tests/io/graph.csv', delimiter=',', comments='#', nodetype=str, data=(('weight', float),('error',float)))

# def test_cycle_closure(pG, capsys):


# def test_compute_weighted_avg_paths():

# def test_compute_average_paths():

# def test_remove_compound_from_graph():

# def test_format_free_energies():
