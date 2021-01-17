import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit.circuit import Parameter

from aqcm_circuits import universal_qcm, universal_qcm_SWAP, phase_covariant_qcm, economical_qcm


def prepare_qubit(qc, qubit, theta, phi):
    """
    Prepares the state |psi>=cos(theta/2)|0>+e^(i phi)sin(theta/2)|1>
    """
    ry = qc.ry(theta, qubit)
    qc.rz(phi, qubit)


def rotated_measurement(qc, qubit, cbit, theta, phi):
    """
    Performs a measurement on qreg qubit (result in creg) in the {|psi>, |psi_perp>} basis,
    where |psi>=cos(theta/2)|0>+e^(i phi)sin(theta/2)|1>
    """

    qc.rz(-phi, qubit)
    qc.ry(-theta, qubit)

    qc.measure(qubit, cbit)


def prepare_qubit_equator(qc, qubit, theta):
    """
    Prepares the state |psi>=cos(theta/2)|0>+sin(theta/2)|1>
    """
    qc.ry(theta, qubit)


def rotated_measurement_equator(qc, qubit, cbit, theta):
    """
    Performs a measurement on qreg qubit (result in creg) in the {|psi>, |psi_perp>} basis,
    where |psi>=cos(theta/2)|0>+sin(theta/2)|1>
    """
    qc.ry(-theta, qubit)

    qc.measure(qubit, cbit)


def sphere_points(num_pts):
    """
    Create num_pts evenly distributed points on a sphere as explained at
    https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    """
    indices = np.arange(0, num_pts, dtype=float) + 0.5

    theta = np.arccos(1 - 2 * indices / num_pts)
    phi = (np.pi * (1 + 5 ** 0.5) * indices) % (2 * np.pi)
    coords = []

    for theta_elem, phi_elem in zip(theta, phi):
        coords.append((theta_elem, phi_elem))

    return np.array(coords)


def equator_xz_points(num_pts):
    """
    Generate num_pts evenly spaced points on the xz equator. It returns their angular coordinates.
    """

    coords = []
    for i in range(num_pts):
        coords.append((2 * np.pi / num_pts * i, 0))

    return np.array(coords)

def prepare_circuits(qcm, input_states, num_points):
    """
    Prepare the circuits for a given qcm, set of input states and number of points (always 4 for the bb84 states)
    It returns the coordinates of the points, the circuits not transpiled and the qubits of the copies.
    qcm can be either 'uqcm', 'uqcm_swap', 'pcqcm', 'epcqcm'.
    input_states can be either 'sphere', 'equator' or 'bb84'.
    bb84 works just like equator, but the number of points is fixed at 4
    """

    # Prepare the target points coordinates
    if input_states == 'sphere':
        target_points=sphere_points(num_points)
    elif input_states == 'equator':
        target_points=equator_xz_points(num_points)
    elif input_states == 'bb84':
        target_points=equator_xz_points(4)
    else:
        print("Error: input states set not recognized.")
        exit()

    # Prepare the parameters for the circuits
    theta_param, phi_param = Parameter('theta_param'), Parameter('phi_param')

    # Define the number of qubits and create the circuits
    if qcm=='uqcm' or qcm=='uqcm_swap' or qcm=='pcqcm':
        nqubits = 3
    elif qcm=='epcqcm':
        nqubits = 2
    else:
        print("Error: quantum cloning machine not recognized.")
        exit()
    qreg = QuantumRegister(nqubits, 'q')
    creg = ClassicalRegister(nqubits, 'c')
    circuit = QuantumCircuit(qreg, creg)

    # Add the preparation to the circuit
    if input_states == 'sphere':
        prepare_qubit(circuit, qreg[0], theta_param, phi_param)
    elif input_states == 'equator' or input_states=='bb84':
        prepare_qubit_equator(circuit, qreg[0], theta_param)
    else:
        print("Error: input states set not recognized.")
        exit()

    # Add the cloning machine circuit block
    if qcm=='uqcm':
        qubit_copy1, qubit_copy2 = universal_qcm(circuit, qreg[0], qreg[1], qreg[2])
    elif qcm=='uqcm_swap':
        qubit_copy1, qubit_copy2 = universal_qcm_SWAP(circuit, qreg[0], qreg[1], qreg[2])
    elif qcm=='pcqcm':
        qubit_copy1, qubit_copy2 = phase_covariant_qcm(circuit, qreg[0], qreg[1], qreg[2])
    elif qcm=='epcqcm':
        qubit_copy1, qubit_copy2 = economical_qcm(circuit, qreg[0], qreg[1])
    else:
        print("Error: quantum cloning machine not recognized.")
        exit()
    index_copy1 = qubit_copy1.index
    index_copy2 = qubit_copy2.index

    # Add the fidelity measurement
    if input_states == 'sphere':
        rotated_measurement(circuit, qreg[index_copy1], creg[index_copy1], theta_param, phi_param)
        rotated_measurement(circuit, qreg[index_copy2], creg[index_copy2], theta_param, phi_param)
        circuits = [circuit.bind_parameters({theta_param: points[0], phi_param: points[1]}) for points in target_points]
    elif input_states == 'equator' or input_states == 'bb84':
        rotated_measurement_equator(circuit, qreg[index_copy1], creg[index_copy1], theta_param)
        rotated_measurement_equator(circuit, qreg[index_copy2], creg[index_copy2], theta_param)
        circuits = [circuit.bind_parameters({theta_param: points[0]}) for points in target_points]
    else:
        print("Error: input states set not recognized.")
        exit()

    # Return the coordinates of the points, the circuits and the qubits of the copies
    return target_points, circuits, qubit_copy1, qubit_copy2


def prepare_circuits_and_transpile(qcm, input_states, num_points, backend, optimization_level=1, transpiler_seed=None, initial_layout=None):
    """
    Transpiling a circuit before binding the parameters might be faster.
    Prepare the circuits for a given qcm, set of input states and number of points (always 4 for the bb84 states)
    It returns the coordinates of the points, the circuits transpiled and the indices of the transpiled copies.
    qcm can be either 'uqcm', 'uqcm_swap', 'pcqcm', 'epcqcm'.
    input_states can be either 'sphere', 'equator' or 'bb84'.
    bb84 works just like equator, but the number of points is fixed at 4
    """

    # Prepare the target points coordinates
    if input_states == 'sphere':
        target_points=sphere_points(num_points)
    elif input_states == 'equator':
        target_points=equator_xz_points(num_points)
    elif input_states == 'bb84':
        target_points=equator_xz_points(4)
    else:
        print("Error: input states set not recognized.")
        exit()

    # Prepare the parameters for the circuits
    theta_param, phi_param = Parameter('theta_param'), Parameter('phi_param')

    # Define the number of qubits and create the circuits
    if qcm=='uqcm' or qcm=='uqcm_swap' or qcm=='pcqcm':
        nqubits = 3
    elif qcm=='epcqcm':
        nqubits = 2
    else:
        print("Error: quantum cloning machine not recognized.")
        exit()
    qreg = QuantumRegister(nqubits, 'q')
    creg = ClassicalRegister(nqubits, 'c')
    circuit = QuantumCircuit(qreg, creg)

    # Add the preparation to the circuit
    if input_states == 'sphere':
        prepare_qubit(circuit, qreg[0], theta_param, phi_param)
    elif input_states == 'equator' or input_states=='bb84':
        prepare_qubit_equator(circuit, qreg[0], theta_param)
    else:
        print("Error: input states set not recognized.")
        exit()

    # Add the cloning machine circuit block
    if qcm=='uqcm':
        qubit_copy1, qubit_copy2 = universal_qcm(circuit, qreg[0], qreg[1], qreg[2])
    elif qcm=='uqcm_swap':
        qubit_copy1, qubit_copy2 = universal_qcm_SWAP(circuit, qreg[0], qreg[1], qreg[2])
    elif qcm=='pcqcm':
        qubit_copy1, qubit_copy2 = phase_covariant_qcm(circuit, qreg[0], qreg[1], qreg[2])
    elif qcm=='epcqcm':
        qubit_copy1, qubit_copy2 = economical_qcm(circuit, qreg[0], qreg[1])
    else:
        print("Error: quantum cloning machine not recognized.")
        exit()
    index_copy1 = qubit_copy1.index
    index_copy2 = qubit_copy2.index

    # Add the fidelity measurement
    if input_states == 'sphere':
        rotated_measurement(circuit, qreg[index_copy1], creg[index_copy1], theta_param, phi_param)
        rotated_measurement(circuit, qreg[index_copy2], creg[index_copy2], theta_param, phi_param)
        circuit_transpiled = transpile(circuit, backend=backend, optimization_level=optimization_level, seed_transpiler=transpiler_seed, initial_layout=initial_layout)
        circuits_transpiled = [circuit_transpiled.bind_parameters({theta_param: points[0], phi_param: points[1]}) for points in target_points]
    elif input_states == 'equator' or input_states == 'bb84':
        rotated_measurement_equator(circuit, qreg[index_copy1], creg[index_copy1], theta_param)
        rotated_measurement_equator(circuit, qreg[index_copy2], creg[index_copy2], theta_param)
        circuit_transpiled = transpile(circuit, backend=backend, optimization_level=optimization_level, seed_transpiler=transpiler_seed, initial_layout=initial_layout)
        circuits_transpiled = [circuit_transpiled.bind_parameters({theta_param: points[0]}) for points in target_points]
    else:
        print("Error: input states set not recognized.")
        exit()

    # Get the qubit index of the copies in the transpiled circuits
    # Need to distinguish between simulator and non-simulator, otherwise error
    if backend.configuration().simulator:
        index_copy1_transpiled = index_copy1
        index_copy2_transpiled = index_copy2
    else:
        index_copy1_transpiled = circuits_transpiled[0].__dict__['_layout'][qubit_copy1]
        index_copy2_transpiled = circuits_transpiled[0].__dict__['_layout'][qubit_copy2]

    # Return the coordinates of the points, the circuits and the indices of the transpiled copies
    return target_points, circuits_transpiled, index_copy1, index_copy2, index_copy1_transpiled, index_copy2_transpiled
