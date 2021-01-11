from qiskit import IBMQ
from quantuminspire.api import QuantumInspireAPI
from analyze_data import *
from execution import *
import numpy as np
import sys

def recover_data(path, target_points):
    """
    Sometimes something happens and a long job is interrupted.
    This function allows to continue the job where it stoped.
    It returns the remaining target points, as well as the intermediate results with and without correction.
    """

    #Specify the files that have partial data
    job = 'results_starmon5.txt'
    job_c = 'results_corrected_starmon5.txt'

    #Recover data
    num_lines = sum(1 for line in open(job))
    results = np.loadtxt(path+job)[:,2:]
    results_c = np.loadtxt(path+job_c)[:,2:]

    #Make output data
    target_points = target_points[num_lines:]
    recovered = [results.tolist(), results_c.tolist()]

    return target_points, recovered

def main():
    """
    Main execution of the code.
    """

    num_pts = 1000
    path = './starmon_bb84/'
    backend_identifier, backend = get_backend(sys.argv[1])

    target_points = bb84_points()
    only_equator = True

    #If you want to recover an interrupted job, set this variable to True
    recover_job = False

    if recover_job == True:
        target_points, recovered = recover_data(path, target_points)
    else:
        recovered = None

    if backend_identifier not in ['starmon5', 'qi_simulator']:
        IBM_execution(only_equator, target_points, backend, backend_identifier, num_pts, path, recovered)

    elif backend_identifier in ['starmon5', 'qi_simulator']:
        Quantum_Inspire_Execution(path, backend, backend_identifier, 1024, target_points, [1,2,3], [1,2], recovered)

    else:
        print('Invalid backend.')

if __name__ == '__main__':
    main()
