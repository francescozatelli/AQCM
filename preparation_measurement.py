import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, assemble, transpile

TOKEN_QI = 'YOUR_TOKEN'
TOKEN_IBMQ = 'YOUR_TOKEN '


def prepare_qubit(qc, qubit, theta, phi):
    """
    Prepares the state |psi>=cos(theta/2)|0>+e^(i phi)sin(theta/2)|1>
    """
    qc.ry(theta, qubit)
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


def prepare_qubit_bb84(qc, qubit, index):
    """
    Prepare a bb84 state with the following convention
    index=0 -> |0>
    index=1 -> |1>
    index=2 -> |+>
    index=3 -> |->
    """
    if index == 1:
        qc.x(qubit)
    elif index == 2:
        qc.sxdg(qubit)
    elif index == 3:
        qc.sx(qubit)


def rotated_measurement_bb84(qc, qubit, cbit, index):
    """
    Performs a bb84 state measurement with the following convention
    index=0 -> |0>
    index=1 -> |1>
    index=2 -> |+>
    index=3 -> |->
    """
    if index == 1:
        qc.x(qubit)
    elif index == 2:
        qc.sx(qubit)
    elif index == 3:
        qc.sxdg(qubit)
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


def equator_xy_points(num_pts):
    """
    Generate num_pts evenly spaced points on the xy equator. It returns their angular coordinates.
    """

    coords = []
    for i in range(num_pts):
        coords.append((np.pi / 2, 2 * np.pi / num_pts * i))

    return np.array(coords)


def bb84_points():
    """
    Generate the coordinates of the bb84 states
    """
    coords = [(0, 0), (np.pi / 2, 0), (np.pi, 0), (3 * np.pi / 2, 0)]

    return np.array(coords)


def calculate_probabilities(job, index_copy1, index_copy2, nshots):
    """
    Calculate the marginal probabilities of getting 0 on each of the two copies qubits
    """
    # Format of the results:
    # {'000': number_of_000, '010': number_of_010, '100': number_of_100, '110': number_of_110}
    marg_prob0_copy1 = 0
    marg_prob0_copy2 = 0
    for outcome in job:
        # If the outcome of the index_copy1 qubit was zero, then add the probability
        if outcome[-1 - index_copy1] == '0':
            marg_prob0_copy1 += job[outcome] / nshots
        # If the outcome of the index_copy2 qubit was zero, then add the probability
        if outcome[-1 - index_copy2] == '0':
            marg_prob0_copy2 += job[outcome] / nshots
    return marg_prob0_copy1, marg_prob0_copy2


def get_counts_from_jobs(jobs):
    """
    Extract the counts from the input jobs
    """
    counts_array = []
    for job in jobs:
        # get_counts() returns a list of Counts (if more experiment were run in the job)
        # or a single object Counts (if only one experiment was run in the job).
        # In order to concatenate (i.e. use the + operator) we have to make sure that job.result().get_counts() is a list
        if isinstance(job.result().get_counts(), list):
            counts_array = counts_array + job.result().get_counts()
        else:
            counts_array = counts_array + [job.result().get_counts()]
    #print(counts_array)
    return counts_array


def analyze_data(job_to_analyze, index_copy1, index_copy2, nshots):
    """
    Returns an array of marginal probabilities for the experiments run in job_to_analyze
    """
    # Get the counts in each job and merge everything into only one array
    results_probabilities = []
    counts_array = get_counts_from_jobs([job_to_analyze])
    # Calculate the marginal probabilities from the counts histograms.
    for job in counts_array:
        marg_prob0_copy1, marg_prob0_copy2 = calculate_probabilities(job, index_copy1, index_copy2,
                                                                     nshots)
        marg_prob1_copy1 = 1 - marg_prob0_copy1
        marg_prob1_copy2 = 1 - marg_prob0_copy2
        results_probabilities.append(
            [marg_prob0_copy1, marg_prob1_copy1, marg_prob0_copy2, marg_prob1_copy2])

    return results_probabilities




