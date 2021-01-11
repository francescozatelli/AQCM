import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, assemble, transpile
from qiskit import IBMQ
from quantuminspire.qiskit import QI
from quantuminspire.api import QuantumInspireAPI
from preparation_measurement import *

def make_readout_circuits():
    """Single qubit calibration run independently for all the qubits in the circuit"""
    
    circuits = []
    n = 2

    for index in range(n):
        qc_0 = QuantumCircuit(n,1)
        qc_1 = QuantumCircuit(n,1)
        qc_1.x(index)
        qc_0.measure(index,0)
        qc_1.measure(index,0)
        circuits.append([qc_0,qc_1])
        
    return circuits

def readout_analyze_data(job_to_analyze, nshots):
    """
    Returns an array of marginal probabilities for the experiments run in job_to_analyze
    """
    # Get the counts in each job and merge everything into only one array
    correcting_parameters = []
    counts_array = get_counts_from_jobs([job_to_analyze])
    reorder_count_array =[counts_array[2*i:(i)*2+2] for i in range(0,len(counts_array)//2)]
    # Calculate the marginal probabilities from the counts histograms.
    for job in reorder_count_array:
        correcting_parameters.append(readout_correction(job, nshots))
    #print(correcting_parameters)

    return correcting_parameters

def get_parameters(p0, p1):
    """
    Calculate coefficients beta0 and beta1 for classical readout correction
    """
    beta_1 = 0.5*(p0+p1)
    beta_2 = 0.5*(p0-p1)

    return beta_1, beta_2

def readout_correction(job, nshots):
    """
    Returns an array of marginal probabilities for the experiments run in job_to_analyze
    """
    
    p0_0 = job[0]['0']/nshots
    mz_0 = p0_0 if p0_0 == 1.0 else p0_0 - job[0]['1']/nshots


    p1_1 = job[1]['1']/nshots
    mz_1 = -p1_1 if p1_1 == 1.0 else -p1_1 + job[1]['0']/nshots


    #print(job.result().get_counts())
    return get_parameters(mz_0,mz_1)


def correct_copies(b,p):
    """
    The results for the probabilities p_0 and p_1 are corrected using the readout parameters b.
    """
    corrected_results = []
    print(b,p)

    #print(p)

    for i in range(len(b)):
        q0_corr= (b[i][1]-b[i][0] + p[i][0] -  p[i][1])/(2*b[i][1])
        q1_corr= 1-q0_corr
        corrected_results += [q0_corr] + [q1_corr]

    #print(corrected_results)

    return corrected_results

def run_readout_correction(readout_obj,nshots, backend):
    print("Readout calibration...")
    readout_params = readout_analyze_data(backend.run(readout_obj),nshots)
    return readout_params
    
def make_readout_circuits_ibm(backend, layout, copies, only_equator, target_points):

    readout_circuits = make_readout_circuits()
    print(readout_circuits)
    readout_circuits_transpiled = transpile(readout_circuits[0]+readout_circuits[1], backend=backend, optimization_level=3,initial_layout=[layout[copies[0]], layout[copies[1]]])
    return readout_circuits_transpiled

def beta2_circuit(qubit):
    qasm = '''
    version 1.0

    qubits 5
    
    X q[{0}]
    
    Measure_z q[{0}]
    '''
    return qasm.format(qubit)

def beta1_circuit(qubit):
    qasm = '''
    version 1.0

    qubits 5

    Measure_z q[{0}]
    '''
    return qasm.format(qubit)

def run_readout_qi(qubit, backend, N_shots):
    qi = QuantumInspireAPI()
    result1 = qi.execute_qasm(beta1_circuit(qubit), number_of_shots=N_shots, backend_type=backend)
    result2 = qi.execute_qasm(beta2_circuit(qubit), number_of_shots=N_shots, backend_type=backend)
    print(result1.get('histogram',{}))
    print(result2.get('histogram',{}))
    measurement_qubit = 2**qubit
    return get_parameters(result1.get('histogram',{})['0'], result2.get('histogram',{})[str(measurement_qubit)])

def readout_correction_qi(index, layout_measurement, limit, backend, N_shots):
    readout_params = []

    for qubit in layout_measurement:
        p0, p1 = run_readout_qi(qubit, backend, N_shots)
        readout_params.append([p1,p0])

    return readout_params