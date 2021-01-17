import numpy as np

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
    Get the histograms counts for each experiment
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
    return counts_array


def analyze_data(job_to_analyze, index_copy1, index_copy2, nshots):
    """
    Calculate the marginal probabilities from a single experiment
    """
    # Get the counts in each job and merge everything into only one array
    results_probabilities = []
    counts_array = get_counts_from_jobs([job_to_analyze])
    # Calculate the marginal probabilities from the counts histograms.
    for job in counts_array:
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

        marg_prob1_copy1 = 1 - marg_prob0_copy1
        marg_prob1_copy2 = 1 - marg_prob0_copy2
        results_probabilities.append(
            [marg_prob0_copy1, marg_prob1_copy1, marg_prob0_copy2, marg_prob1_copy2])

    return results_probabilities


def print_single_result(path, backend_identifier, theta, phi, results_probabilities, num_pts, label='', print_on_screen=False):
    """
    Print results both on terminal and on a file
    """
    with open(path + 'results_' + backend_identifier + '_' + label + '.txt', 'a') as file:
        file.write(
            str(theta) + "\t" + str(phi) + "\t" + str(
                results_probabilities[-1][0]) + "\t" + str(
                results_probabilities[-1][1]) + "\t" +
            str(results_probabilities[-1][2]) + "\t" + str(results_probabilities[-1][3]) + "\n")

    if print_on_screen:
        print(len(results_probabilities), "/", num_pts, label)
        print("Coordinates theta,phi: ", [theta, phi])
        print("Fidelities copy1,copy2: ", [results_probabilities[-1][0], results_probabilities[-1][2]])
        print("Average fidelity copy1, copy2: ", [np.average(np.array(results_probabilities)[:, 0]),
                                                  np.average(np.array(results_probabilities)[:, 2])])


def write_average_fidelity(path, backend_identifier, results_probabilities, label=''):
    """
    Save average fidelity in a file
    """
    with open(path + 'average_fidelities_' + label + '.txt', 'a') as file:
        file.write(
            backend_identifier + "\t" + str(np.average(np.array(results_probabilities)[:, 0])) + "\t" + str(
                np.std(np.array(results_probabilities)[:, 0])) + "\t" + str(
                np.average(np.array(results_probabilities)[:, 2])) + "\t" +
            str(np.std(np.array(results_probabilities)[:, 2])) + "\n")



def save_result(theta, phi, result_probabilities, result_probabilities_corrected, backend_identifier, path):
    # Save non corrected results
    with open(path + 'results_' + backend_identifier+ '.txt', 'a') as file:
        file.write(
            str(theta) + "\t" + str(phi) + "\t" + str(
                result_probabilities[0]) + "\t" + str(
                result_probabilities[1]) + "\t" +
            str(result_probabilities[2]) + "\t" + str(result_probabilities[3]) + "\n")

    # Save corrected results
    with open(path + 'results_' + backend_identifier + '_corrected.txt', 'a') as file:
        file.write(
            str(theta) + "\t" + str(phi) + "\t" + str(
                result_probabilities_corrected[0]) + "\t" + str(
                result_probabilities_corrected[1]) + "\t" +
            str(result_probabilities_corrected[2]) + "\t" + str(result_probabilities_corrected[3]) + "\n")