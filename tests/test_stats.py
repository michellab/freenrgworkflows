import pytest
import numpy as np
from networkanalysis.stats import *

@pytest.fixture
def stats():
    return freeEnergyStats()

def test_calculate_r2(stats):
    a = np.array([2,4,5,6,4,6,8.9,7])
    b = np.array([5,4,6,7,8,7,7.8,6])
    r2 = 0.2931274858030396 
    r = 0.5414124913622141
    test_r2, test_r = stats._calculate_r2(a,b)
    assert(test_r2 == r2)
    assert(test_r == r)


def test_calculate_tau(stats):
    a = np.array([2,4,5,6,4,6,8.9,7])
    b = np.array([5,4,6,7,8,7,7.8,6])
    tau = 0.3461538461538462 
    test_tau= stats._calculate_tau(a,b)
    assert(test_tau == tau)

def test_calculate_mue(stats):
    a = np.array([2,4,5,6,4,6,8.9,7])
    b = np.array([5,4,6,7,8,7,7.8,6])
    mue = 1.5125000000000002
    test_mue= stats._calculate_mue(a,b)
    assert(test_mue == mue)

