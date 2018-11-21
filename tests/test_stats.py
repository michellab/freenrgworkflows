import pytest
import numpy as np
import warnings
from networkanalysis.stats import *


@pytest.fixture
def stats():
    return freeEnergyStats()


def test_calculate_r2(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    r2 = 0.2931274858030396
    r = 0.5414124913622141
    test_r2, test_r = stats._calculate_r2(a, b)
    assert (test_r2 == r2)
    assert (test_r == r)


def test_calculate_tau(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    tau = 0.3461538461538462
    test_tau = stats._calculate_tau(a, b)
    assert (test_tau == tau)


def test_calculate_mue(stats):
    a = np.array([2, 4, 5, 6, 4, 6, 8.9, 7])
    b = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    mue = 1.5125000000000002
    test_mue = stats._calculate_mue(a, b)
    assert (test_mue == mue)


@pytest.mark.parametrize('boundaries', [(-72), (4)])
def test_confidence_warnings(stats, boundaries):
    warn_string = 'Confidence interval needs to be between 0 and 1, please try something like 0.68 for one sigma confidence'
    with pytest.warns(UserWarning) as warnmessg:
        stats.confidence_interval = boundaries
        # warnings.warn(warn_string, UserWarning)
    assert len(warnmessg) == 1
    assert warnmessg[0].message.args[0] == warn_string


def test_confidence(stats):
    data = np.array([5, 4, 6, 7, 8, 7, 7.8, 6])
    assert (stats._confidence(data) == [6.0, 7.8])


@pytest.mark.parametrize('repeat', [(100), (20)])
def test_repeats(stats, repeat):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=repeat)
    assert (len(stats._R) == repeat)
    assert (len(stats._mue) == repeat)
    assert (len(stats._tau) == repeat)
    assert (len(stats._R2) == repeat)


def test_array_conversion(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=1)
    data_comp = np.array(stats.data_comp)
    data_exp = np.array(stats.data_exp)
    assert (np.all(data_comp[:, 0]) == np.all(data_exp))


def test_statistics(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=1)
    comp_array = np.array(stats.data_comp)
    r = stats._calculate_r2(comp_array[:, 0], stats.data_exp)
    tau = stats._calculate_tau(comp_array[:, 0], stats.data_exp)
    mue = stats._calculate_mue(comp_array[:, 0], stats.data_exp)
    assert (r == (1.0, 1.0))
    assert (tau == 1.0)
    assert (mue == 0.0)


def test_calculate_predictive_index(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    error_string = 'Calculating predictive index not impletmented yet.'
    with pytest.raises(NotImplementedError, message=error_string):
        stats._calculate_predictive_index(exp_dat, comp)


def test_properties(stats):
    exp_dat = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    comp = [{'a': 0.5, 'error': 0.02}, {'b': 1.7, 'error': 0.2}, {'f': -0.9, 'error': 0.06}]
    stats.generate_statistics(comp, exp_dat, repeats=100)
    assert (np.mean(stats._R) == stats.R)
    assert (np.mean(stats._R2) == stats.R2)
    assert (np.mean(stats._tau) == stats.tau)
    assert (np.mean(stats._mue) == stats.mue)


def test_property_errors(stats):
    stats._R = np.loadtxt('tests/io/R.dat')
    stats.confidence_interval = 0.68
    err1 = stats.R_confidence
    stats.confidence_interval = 0.95
    err2 = stats.R_confidence
    assert (err1 != err2)


@pytest.mark.parametrize('boundaries', [(0.95)])
def test_property_errors(stats, boundaries):
    stats._R = np.loadtxt('tests/io/R.dat')
    stats._tau = np.loadtxt('tests/io/tau.dat')
    stats._mue = np.loadtxt('tests/io/mue.dat')
    stats._R2 = stats._R ** 2
    stats.confidence_interval = boundaries
    assert (stats.R_confidence[1:] == [0.8682078950591153, 0.9338724006518556])
    assert (stats.tau_confidence[1:] == [1.0, 1.0])
    assert (stats.mue_confidence[1:] == [0.3562336854267332, 0.5309578343443963])
    assert (stats.R2_confidence[1:] == [0.7537849490429798, 0.8721176606992599])
