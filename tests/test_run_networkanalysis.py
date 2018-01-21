import subprocess
import pytest
import os
import sys
from pathlib import Path
import subprocess

@pytest.fixture
def executable():
    return os.path.join(os.getcwd(),'bin','run_networkanalysis.py')
@pytest.fixture
def graph_file():
    return os.path.join(os.getcwd(),'tests','io', 'graph.csv')

def test_input_test_data_exists(executable,graph_file):
    assert(Path(executable).exists())
    assert(Path(graph_file).exists())

def test_no_file_passed(executable):
    cmd = [sys.executable, executable]
    error_string = b'You must give at least one networkanalysis networkx compatible file'
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(error_string in stderr)

def test_no_target_compound(executable, graph_file):
    cmd = [sys.executable, executable, graph_file]
    error_string = b'No target compound given, using the first compound in the node list:'
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(error_string in stderr)

def test_wrong_target_compound(executable, graph_file):
    cmd = [sys.executable, executable, graph_file ,'--target_compound=FXR20','--network_output=tests/io/test_out.dat']
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    err_string = b'source node FXR20 not in graph'
    assert(p.returncode==1)
    assert(err_string in stderr)

def test_safe_data(executable, graph_file):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'test_out.dat')
    cmd = [sys.executable, executable, graph_file ,'--target_compound=FXR17','--network_output='+filename]
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(p.returncode == 0)
    assert(Path(filename).exists())
    os.remove(Path(filename))

def test_do_not_safe_data(executable, graph_file):
    cmd = [sys.executable, executable, graph_file ,'--target_compound=FXR17']
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(p.returncode == 0)
    assert(b'#FREE ENERGIES ARE:' in stdout)

def test_statistics(executable, graph_file):
    filename = os.path.join(os.getcwd(), 'tests', 'io', 'ic50_exp.dat')
    cmd = [sys.executable, executable, graph_file ,'--target_compound=FXR17','--experiments='+filename]
    output = b'R and std'
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(p.returncode == 0)
    assert(output in stdout)

    