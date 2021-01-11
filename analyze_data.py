import numpy as np
from preparation_measurement import *
from readout_correction import *
from qiskit import IBMQ
from quantuminspire.credentials import enable_account
from quantuminspire.qiskit import QI
from quantuminspire.api import QuantumInspireAPI

list_backend_identifiers_ibm = ['simulator',
                                'ibmqx2',
                                'melbourne',
                                'vigo',
                                'ourense',
                                'valencia',
                                'athens',
                                'santiago']

list_backend_identifiers_qi = ['starmon5', 'qi_simulator']    

def get_backend(inputt):
    inputt = str(inputt)

    if inputt in list_backend_identifiers_ibm:

        IBMQ.save_account('73b1911cec7f036d146d1900a34a4b29b65288ca26463d7baae30451d376042d501673bd7902118811502a2a02b89126af12c64c14ad5ae094fe36db1a2ed9cd', overwrite=True)
        provider = IBMQ.load_account()
        list_backends = [provider.backends.ibmq_qasm_simulator,
                     provider.backends.ibmqx2,
                     provider.backends.ibmq_16_melbourne,
                     provider.backends.ibmq_vigo,
                     provider.backends.ibmq_ourense,
                     provider.backends.ibmq_valencia,
                     provider.backends.ibmq_athens,
                     provider.backends.ibmq_santiago]

        index = list_backend_identifiers_ibm.index(inputt)
        return inputt, list_backends[index]

    elif inputt in list_backend_identifiers_qi:

        enable_account('55e8b5df17d6f0cbdcc7fec27245b744f0484303')
        QI.set_authentication()
        qi = QuantumInspireAPI()

        list_backends = [qi.get_backend_type_by_name('Starmon-5'),qi.get_backend_type_by_name('QX Single-Node Simulator')]
        index = list_backend_identifiers_qi.index(inputt)

        return inputt, list_backends[index]
    
    else:
        print('Invalid Backend!')

def write_readout_parameters(readout, backend_identifier, path):
    """
    Write results of readout parameters in a file
    """

    with open(path + 'results_readout_' + backend_identifier + '.txt', 'a') as file:

        for i in range(len(readout)):
            for j in range(len(readout[i])):
                file.write(str(readout[i][j])+'\t')
        file.write('\n')


def write_copying_results(coords, results, label, backend_identifier, path):
    """
    Write results of qubit measurements in a file
    """

    with open(path + 'results_' + label + backend_identifier + '.txt', 'a') as file:
        file.write(str(coords[0]) + "\t" + str(coords[1]) + "\t")
        for i in range(len(results[-1])):
            file.write(
                str(results[-1][i]) + "\t")
        file.write('\n')

def write_average_fidelities(results_list, backend_identifier, path):
    """
    Write fidelities for corrected and original results.
    """

    labels = ['','corrected']
    for i in range(len(results_list)):
        with open(path + 'average_fidelities.txt', 'a') as file:
            file.write(
                backend_identifier+"_"+str(labels[i])+"\t" + str(np.average(np.array(results_list[i])[:, 0])) + "\t" + str(
                    np.std(np.array(results_list[i])[:, 0])) + "\t" + str(
                    np.average(np.array(results_list[i])[:, 2])) + "\t" +
                str(np.std(np.array(results_list[i])[:, 2])) + "\n")


def save_experiment(batch, results_probabilities, corrected_results, readout_params, backend_identifier,target_points, path):
    """
    Receives measurement results and readout parameters, correct the measurements
    and write resulting data in files
    """
    #print(batch, readout_params)
    # Append the results and print them one by one
    for result in batch:
        print(results_probabilities)
        results_probabilities.append(result)
        coords = (
            target_points[len(results_probabilities) - 1][0],
            target_points[len(results_probabilities) - 1][1])
                
        #Readout calibration
        #print(results_probabilities)
        reordered_results = [results_probabilities[-1][i:i+2] for i in range(0,len(results_probabilities[-1])-1,2)]
        corrected_results.append(correct_copies(readout_params, reordered_results))

        #Write corrected results
        write_copying_results(coords, corrected_results,'corrected_',backend_identifier, path)
        write_copying_results(coords, results_probabilities,'',backend_identifier, path)
                
        #print(len(results_probabilities), "/", num_pts)
    
    return results_probabilities, corrected_results


def get_layout_measurement_qubits(backend_identifier, circuits_transpiled, nqubits=3):
    #Find the positions of the physical qubits transpiled for the circuit
    layout = []

    if backend_identifier != 'simulator':
        for element in circuits_transpiled[0].__dict__['_layout'].get_physical_bits():
            layout.append(element)

    else:
        layout = range(nqubits)

    return layout

def get_marginal_probabilities_qi(hist):
    #test both the copies

    if hist['0'] == 1.0:
        marg_prob0_qubit0=hist['0']
        marg_prob0_qubit1=hist['0']
    
    else:
        marg_prob0_qubit0=hist['0']+hist['4'] #calculate the marginal distribution, p=p(00000)+p(00010)
        marg_prob0_qubit1=hist['0']+hist['2'] #calculate the marginal distribution, p=p(00000)+p(00001)

    marg_prob1_qubit0=1-marg_prob0_qubit0
    marg_prob1_qubit1=1-marg_prob0_qubit1

    data = [marg_prob0_qubit0, marg_prob1_qubit0, marg_prob0_qubit1, marg_prob1_qubit1]
    print(data)
    return data

def get_correct_results_id(layout):
    pass
