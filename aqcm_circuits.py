import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, assemble, transpile
from qiskit.circuit import Parameter
from preparation_measurement import *
from readout_correction import *

def universal_qcm(qc, input, ancilla1, ancilla2):
    """
    Basic circuit for universal cloning machine. The copies will be on the input qubit and on the first ancilla qubit.
    The qubits of the two copies are returned.
    No SWAP gate is used here (all the qubits must be connected).
    Ideal fidelity F=5/6=0.8333...
    """
    # Rotation angles used to prepare the ancillae
    theta1 = np.arccos(1 / np.sqrt(5))
    theta2 = np.arccos(np.sqrt(5) / 3)
    theta3 = np.arccos(2 / np.sqrt(5))

    # Preparation of the ancillae
    qc.ry(theta1, ancilla1)
    qc.cx(ancilla1, ancilla2)
    qc.ry(theta2, ancilla2)
    qc.cx(ancilla2, ancilla1)
    qc.ry(theta3, ancilla1)

    # Actual copying
    qc.cx(input, ancilla1)
    qc.cx(input, ancilla2)
    qc.cx(ancilla1, input)
    qc.cx(ancilla2, input)

    return input, ancilla1


def universal_qcm_SWAP(qc, input, ancilla1, ancilla2):
    """
    Basic circuit for universal cloning machine. The copies will be on the input qubit and on the first ancilla qubit.
    With this circuit only the first ancilla has to be connected to both the other qubits.
    A swap gate is then used in order to perform two qubits gates between the input state and both the first ancilla.
    The qubits of the two copies are returned.
    This circuit can be used even if not all the qubits are mutually connected.
    Ideal fidelity F=5/6=0.8333...
    """

    # Rotation angles used to prepare the ancillae
    theta1 = np.arccos(1 / np.sqrt(5))
    theta2 = np.arccos(np.sqrt(5) / 3)
    theta3 = np.arccos(2 / np.sqrt(5))

    # Preparation of the ancillae
    qc.ry(theta1, ancilla1)
    qc.cx(ancilla1, ancilla2)
    qc.ry(theta2, ancilla2)
    qc.cx(ancilla2, ancilla1)
    qc.ry(theta3, ancilla1)

    # Swap the input qubit and the ancilla1 qubit
    qc.swap(input, ancilla1)

    # Actual copying
    qc.cx(ancilla1, input)
    qc.cx(ancilla1, ancilla2)
    qc.cx(input, ancilla1)
    qc.cx(ancilla2, ancilla1)

    return input, ancilla1

def parameterized_QACM(loc, layout):
    θ,ϕ=loc
    q1,q2,q3 = layout
    qasm = '''
    version 1.0

    qubits 5

    # initialize the state
    Ry q[{2}], {0}
    Rz q[{2}], {1}

    #preparation
    Ry q[{3}], 1.107149
    #rewrite CNOT q[0],q[4] and CNOT q[4],q[0] usign nearest neighbors
    CNOT q[{3}],q[{4}]

    Ry q[{4}], 0.729728
    CNOT q[{4}],q[2]

    Ry q[{3}], 0.463648
    SWAP q[{2}],q[{3}]

    #copying
    CNOT q[{3}], q[{2}]
    CNOT q[{3}], q[{4}]
    CNOT q[{2}], q[{3}]
    CNOT q[{4}], q[{3}]

    #Rotate back and measure
    Rz q[{3}], -{1}
    Ry q[{3}], -{0}
    Rz q[{2}], -{1}
    Ry q[{2}], -{0}
    Measure_z q[{3}]
    Measure_z q[{2}]
    '''
    return qasm.format(θ,ϕ,q1,q2,q3)

def phase_covariant_qcm(qc, input, ancilla1, ancilla2):
    """
    Phase covariant cloning machine. Optimally clones the states on the xz equator.
    The copies are made on the two ancilla qubits, which are returned. Only the input qubit needs to be connected to the other two.
    Ideal fidelity F=1/2+1/sqrt(8)=0.854...
    """
    # Rotation angle used to prepare the ancillae
    theta = np.pi / 4

    # Preparation of the ancillae
    qc.ry(theta, ancilla1)
    qc.ry(theta, ancilla2)

    # Actual copying
    qc.cx(input, ancilla2)
    qc.cx(input, ancilla1)
    qc.cx(ancilla2, input)
    qc.cx(ancilla1, input)

    return ancilla1, ancilla2


def prepare_circuits_ibm(only_equator, target_points, backend):

    theta_param, phi_param = Parameter('theta_param'), Parameter('phi_param')
    nqubits = 3
    qreg = QuantumRegister(nqubits, 'q')
    creg = ClassicalRegister(nqubits, 'c')
    circuit = QuantumCircuit(qreg, creg)

    if not only_equator:
        prepare_qubit(circuit, qreg[0], theta_param, phi_param)
    else:
        prepare_qubit_equator(circuit, qreg[0], theta_param)

    qubit_copy1, qubit_copy2 = universal_qcm(circuit, qreg[0], qreg[1], qreg[2])
    index_copy1 = qubit_copy1.index
    index_copy2 = qubit_copy2.index

    if not only_equator:
        rotated_measurement(circuit, qreg[index_copy1], creg[index_copy1], theta_param, phi_param)
        rotated_measurement(circuit, qreg[index_copy2], creg[index_copy2], theta_param, phi_param)
    else:
        rotated_measurement_equator(circuit, qreg[index_copy1], creg[index_copy1], theta_param)
        rotated_measurement_equator(circuit, qreg[index_copy2], creg[index_copy2], theta_param)

    # Prepare circuits
    if not only_equator:
        circuits = [circuit.bind_parameters({theta_param: points[0], phi_param: points[1]}) for points in target_points]
    else:
        circuits = [circuit.bind_parameters({theta_param: points[0]}) for points in target_points]

    circuits_transpiled = transpile(circuits, backend=backend, optimization_level=3)

    return circuits_transpiled, index_copy1, index_copy2