import subprocess
import pytest
import os
import sys

def test_no_file_passed():
    executable = os.path.join(os.getcwd(),'bin','run_networkanalysis.py')
    assert(os.path.exists(executable)==True)
    #cmd = 'python '+executable+' '+graph_file
    cmd = [sys.executable, executable]
    error_string = b'You must give at least one networkanalysis networkx compatible file'
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(error_string in stderr)

def test_no_target_compound():
    executable = os.path.join(os.getcwd(),'bin','run_networkanalysis.py')
    graph_file = os.path.join(os.getcwd(),'tests','io', 'graph.csv')
    assert(os.path.exists(executable)==True)
    assert(os.path.exists(graph_file)==True)
    #cmd = 'python '+executable+' '+graph_file
    cmd = [sys.executable, executable, graph_file]
    error_string = b'No target compound given, using the first compound in the node list:'
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    assert(error_string in stderr)

#def test_wrong_target_compound():
#    executable = os.path.join(os.getcwd(),'bin','run_networkanalysis.py')
#    graph_file = os.path.join(os.getcwd(),'tests','io', 'graph.csv')
#    assert(os.path.exists(executable)==True)
#    assert(os.path.exists(graph_file)==True)
    #cmd = 'python '+executable+' '+graph_file
#    cmd = [sys.executable, executable, graph_file, " --target_compound FXR20" ]
#    error_string = b'node FXR20 not in graph'
#    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
#    stdout, stderr = p.communicate()
#    print (stderr)
    #assert(error_string in stderr)

#def test_safe_data():
#    executable = os.path.join(os.getcwd(),'bin','run_networkanalysis.py')
#    graph_file = os.path.join(os.getcwd(),'tests','io', 'graph.csv')
#    assert(os.path.exists(executable)==True)
#    assert(os.path.exists(graph_file)==True)
#    cmd = [sys.executable, executable, graph_file]
#    error_string = b'No target compound given, using the first compound in the node list.'
#    p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
#    stdout, stderr = p.communicate()

