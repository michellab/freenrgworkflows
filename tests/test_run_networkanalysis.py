import subprocess
import pytest
import os
import sys

def test_target_compound():
    executable = os.path.join(os.getcwd(),'bin','run_networkanalysis.py')
    graph_file = os.path.join(os.getcwd(),'tests','io', 'graph.csv')
    assert(os.path.exists(executable)==True)
    assert(os.path.exists(graph_file)==True)
    #cmd = 'python '+executable+' '+graph_file
    cmd = [sys.executable, executable]
    error_string = 'you must give at least one networkanalysis networkx compatible file'
    #with pytest.raises(IOError, message = error_string):
    #    subprocess.check_output(cmd)