from qiskit import ClassicalRegister, QuantumRegister, assemble, QuantumCircuit, transpile

from data_analysis import get_counts_from_jobs


def generic_readout_circuits():
    """
    Single qubit calibration
    """

    qr_0 = QuantumRegister(1, 'q')
    cr_0 = ClassicalRegister(1)
    qr_1 = QuantumRegister(1, 'q')
    cr_1 = ClassicalRegister(1)

    qc_0 = QuantumCircuit(qr_0, cr_0)
    qc_1 = QuantumCircuit(qr_1, cr_1)

    qc_1.x(qr_1[0])

    qc_0.measure(qr_0[0], cr_0[0])
    qc_1.measure(qr_1[0], cr_1[0])

    return [qc_0, qc_1]


def calculate_readout_parameters_from_probabilities(p0_0, p0_1):
    """
    p0_0 is the probability of outcome 0 in the first readout circuit, p0_1 is the probability of outcome 1 in the second readout circuit
    """

    # Calculate probailities for 1
    p1_0 = 1 - p0_0
    p1_1 = 1 - p0_1

    # Calculate averages
    mA = p0_0 - p1_0
    mB = p0_1 - p1_1

    # Calculate correction parameters
    beta0 = 0.5 * (mA + mB)
    beta1 = 0.5 * (mA - mB)
    return beta0, beta1


def readout_analyze_data(readout_jobs, nshots):
    """
    readout_jobs is the array of the 4 readout jobs [job1, job2, job3, job4]
    """
    counts_array = get_counts_from_jobs(readout_jobs)

    # probability of getting 0 for the first circuit
    p0_0_qubit0 = counts_array[0]['0'] / nshots
    # probability of getting 1 for the second circuit
    p1_1_qubit0 = counts_array[1]['1'] / nshots
    # probability of getting 0 for the second circuit (this is necessary because in the simulator you never get 0 in the second circuit
    p0_1_qubit0 = 1-p1_1_qubit0
    #Same thing
    p0_0_qubit1 = counts_array[2]['0'] / nshots
    p1_1_qubit1 = counts_array[3]['1'] / nshots
    p0_1_qubit1 = 1-p1_1_qubit1

    beta0_qubit0, beta1_qubit0 = calculate_readout_parameters_from_probabilities(p0_0_qubit0, p0_1_qubit0)
    beta0_qubit1, beta1_qubit1 = calculate_readout_parameters_from_probabilities(p0_0_qubit1, p0_1_qubit1)

    return [beta0_qubit0, beta1_qubit0, beta0_qubit1, beta1_qubit1]


def make_readout_circuits(index_physical_qubit, backend, nshots):
    readout_circuits = generic_readout_circuits()
    readout_circuits_transpiled = transpile(readout_circuits, backend=backend, optimization_level=0,
                                            initial_layout=[index_physical_qubit])
    return readout_circuits_transpiled


def readout_calibration(index_physical_qubit, backend, nshots):
    readout_circuits = generic_readout_circuits()
    readout_circuits_transpiled = transpile(readout_circuits, backend=backend, optimization_level=0,
                                            initial_layout=[index_physical_qubit])
    readout_obj = assemble(readout_circuits_transpiled, backend=backend, shots=nshots)
    readout_params = readout_analyze_data(backend.run(readout_obj), nshots)
    return readout_params


def correct_probabilities(result, readout_params):
    beta0_qubit0 = readout_params[0]
    beta1_qubit0 = readout_params[1]
    beta0_qubit1 = readout_params[2]
    beta1_qubit1 = readout_params[3]

    marg_prob0_qubit0 = result[0]
    marg_prob1_qubit0 = result[1]
    marg_prob0_qubit1 = result[2]
    marg_prob1_qubit1 = result[3]

    marg_prob0_qubit0_corr = (beta1_qubit0 - beta0_qubit0 + marg_prob0_qubit0 - marg_prob1_qubit0) / (2 * beta1_qubit0)
    marg_prob1_qubit0_corr = 1 - marg_prob0_qubit0_corr
    marg_prob0_qubit1_corr = (beta1_qubit1 - beta0_qubit1 + marg_prob0_qubit1 - marg_prob1_qubit1) / (2 * beta1_qubit1)
    marg_prob1_qubit1_corr = 1 - marg_prob0_qubit1_corr

    return [marg_prob0_qubit0_corr, marg_prob1_qubit0_corr, marg_prob0_qubit1_corr, marg_prob1_qubit1_corr]


def save_readout_parameters(path, backend_identifier, readout_parameters, index_copy1_transpiled,
                            index_copy2_transpiled):
    with open(path + 'readout_correction_' + backend_identifier + '.txt', 'a') as readout:
        readout.write(
            "Readout correction on qubits: " + str(index_copy1_transpiled) + " and " + str(index_copy2_transpiled) + "\n")
        readout.write(str(readout_parameters[0]) + "\t" + str(readout_parameters[1]) + "\t" + str(
            readout_parameters[2]) + "\t" + str(readout_parameters[3]))
        readout.write('\n')
