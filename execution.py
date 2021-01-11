from preparation_measurement import *
from aqcm_circuits import *
from readout_correction import *
from qiskit import IBMQ
from qiskit.transpiler import Layout
from getpass import getpass
from coreapi.auth import BasicAuthentication
from quantuminspire.api import QuantumInspireAPI
from analyze_data import *
import sys

def IBM_execution(only_equator, target_points, backend, backend_identifier, num_pts, path, recover_job):

    running_jobs = []
    results_probabilities = []
    corrected_results = []
    index = 0

    if recover_job == None:
        results_probabilities = []
        corrected_results = []

    else:
        results_probabilities = recover_job[0]
        corrected_results = recover_job[1]

    max_shots = backend.configuration().max_shots
    max_experiments = backend.configuration().max_experiments

    circuits_transpiled, index_copy1, index_copy2 = prepare_circuits_ibm(only_equator, target_points, backend)
    layout = get_layout_measurement_qubits(backend_identifier, circuits_transpiled)
    print(layout)
    readout_circuits_transpiled = make_readout_circuits_ibm(backend, layout, [index_copy1, index_copy2], only_equator, target_points)
    readout_obj = assemble(readout_circuits_transpiled, backend=backend, shots=max_shots)

    while index * max_experiments < num_pts:

        # Split the transpiled circuits array in an array that contains the maximum number of circuits allowed to run at
        # the same time
        first_circuit_index = index * max_experiments
        last_circuit_index = np.minimum((index + 1) * max_experiments, num_pts)
        max_circuits_transpiled = [circuit for circuit in circuits_transpiled[first_circuit_index:last_circuit_index]]

        #Run a readout before the first round of experiments and calculate data
        readout_params = run_readout_correction(readout_obj, max_shots, backend)
        write_readout_parameters(readout_params, backend_identifier,path)

        # Create and run the maximum number of experiments
        qobj = assemble(max_circuits_transpiled, backend=backend, shots=max_shots)

        #Run the circuits and save the job in an arrray
        running_jobs.append(backend.run(qobj))

        #Save and analyze data from last job
        print(backend_identifier+": Waiting for a job to finish...")
        running_jobs[0].result()  # Wait for the first job to finish
        results_probabilities_batch = analyze_data(running_jobs[0], index_copy1, index_copy2, max_shots)
        results_probabilities, corrected_results = save_experiment(results_probabilities_batch, results_probabilities, corrected_results, readout_params, backend_identifier, target_points, path) 

        #Eliminate the job from the list and get next one
        running_jobs.pop(0)
        index = index + 1

        #If it is the last cycle, we evaluate remaining jobs similarly as before
        if index * max_experiments > num_pts:
            while len(running_jobs) != 0:
                results_probabilities_batch = analyze_data(running_jobs[0], index_copy1, index_copy2, max_shots)
                results_probabilities, corrected_results = save_experiment(results_probabilities_batch, results_probabilities, corrected_results, readout_params, backend_identifier, target_points, path) 
                running_jobs.pop(0)

    #print([results_probabilities, corrected_results])
    write_average_fidelities([results_probabilities, corrected_results],backend_identifier, path)


def Quantum_Inspire_Execution(path, backend, backend_identifier, N_shots, target_points, layout, layout_measurement, recover_job):
    """
    Code that executes a given circuit in Quantum Inspire and writes the data similarly as in IBM_execution.
    
    Input variables:
        path: where files will be stored
        backend, backend_identifier: backend object and label
        N_shots: number of shots
        target_points: points where the circuit will be evaluated
        layout: list with the qubits involved in the circuit
        layour_measurement: list with the qubits that will be measured in the circuit
        recover_job: data regarding a stopped job that is being recovered
    
    """

    qi = QuantumInspireAPI()
    readout_params = []

    if recover_job == None:
        results=[]
        corrected_results=[]

    else:
        results = recover_job[0]
        corrected_results = recover_job[1]

    for index in range(len(target_points)):
        
        if index %75 == 0:
            readout_params_run = readout_correction_qi(index, layout_measurement, 75, backend, N_shots)
            readout_params.append(readout_params_run)
        
        points = target_points[index]
        qasm = parameterized_QACM(points, layout)
        result = qi.execute_qasm(qasm, backend_type=backend, number_of_shots=N_shots)
        batch = get_marginal_probabilities_qi(result.get('histogram',{}))
        
        
        results, corrected_results = save_experiment([batch], results, corrected_results, readout_params_run, backend_identifier, target_points, path) 
    
    #print(results,corrected_results)
    write_average_fidelities([results, corrected_results],backend_identifier, path)

