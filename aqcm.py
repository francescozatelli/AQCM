import sys
from qiskit import assemble
from qiskit import transpile
from qiskit.qobj import QobjHeader
from queue import Queue

from quantuminspire.qiskit import QI

from backend_helper import find_backend_by_nickname
from data_analysis import analyze_data, save_result, write_average_fidelity
from preparation_measurement import prepare_circuits_and_transpile
import numpy as np
from readout_correction import make_readout_circuits, readout_analyze_data, save_readout_parameters, \
    correct_probabilities


def main():
    backend_identifier = sys.argv[
        1]  # This string is used to save all the different files and to find the corresponding backend
    backend = find_backend_by_nickname(backend_identifier)

    max_shots = backend.configuration().max_shots
    max_experiments = backend.configuration().max_experiments
    is_simulator = backend.configuration().simulator
    readout_calibration_step = 75

    """
    Prepare the circuits for a given qcm, set of input states and number of points (always 4 for the bb84 states)
    It returns the coordiantes of the points, the circuits not transpiled and the indices of the copies.
    qcm can be either 'uqcm', 'uqcm_swap', 'pcqcm', 'epcqcm'.
    input_states can be either 'sphere', 'equator' or 'bb84.
    """
    path = "UQCM/IBM/OnlyEquator/"  # Path where the files will be saved (set whatever you want, create folder before using it)
    qcm = 'uqcm_swap'
    input_states = 'equator'
    num_pts = 100

    # It is possible to use a specific layout:
    # layout=[1,2,3]
    # target_points, circuits_transpiled, index_copy1, index_copy2, index_copy1_transpiled, index_copy2_transpiled = prepare_circuits_and_transpile(qcm, input_states, num_pts, backend=backend,optimization_level=2,initial_layout=layout)

    target_points, circuits_transpiled, index_copy1, index_copy2, index_copy1_transpiled, index_copy2_transpiled = prepare_circuits_and_transpile(
        qcm, input_states, num_pts, backend=backend, optimization_level=2)

    # Save the target_points
    np.savetxt(path + "target_points_" + backend_identifier + ".csv", target_points)

    print("Example transpiled circuit", backend_identifier)
    print(circuits_transpiled[0])

    # Prepare readout circuits
    readout_circuits_copy1 = make_readout_circuits(index_copy1_transpiled, backend, max_shots)
    readout_circuits_copy2 = make_readout_circuits(index_copy2_transpiled, backend, max_shots)
    readout_circuits = [readout_circuits_copy1[0], readout_circuits_copy1[1], readout_circuits_copy2[0],
                        readout_circuits_copy2[1]]

    # Prepare qobj queue
    qobj_queue = Queue()
    index = 0
    while index < num_pts:
        if index % readout_calibration_step == 0:
            if backend.provider() == QI:
                # Quantum Inspire does not allow to run more than 3 circuits at a time, so it is easier to just run the readout circuits one by one
                qobj_queue.put(assemble(readout_circuits[0], backend=backend, shots=max_shots,
                                        qobj_header=QobjHeader(is_readout=True, readout_index=0)))
                qobj_queue.put(assemble(readout_circuits[1], backend=backend, shots=max_shots,
                                        qobj_header=QobjHeader(is_readout=True, readout_index=1)))
                qobj_queue.put(assemble(readout_circuits[2], backend=backend, shots=max_shots,
                                        qobj_header=QobjHeader(is_readout=True, readout_index=2)))
                qobj_queue.put(assemble(readout_circuits[3], backend=backend, shots=max_shots,
                                        qobj_header=QobjHeader(is_readout=True, readout_index=3)))
            else:
                # For IBMQ it's better to keep more circuits in a single job, since we can run 5 jobs at a time (each job with 75 circuits)
                qobj_queue.put(assemble(readout_circuits, backend=backend, shots=max_shots,
                                        qobj_header=QobjHeader(is_readout=True)))

        # Put the actual circuits in queue. Add indices and coordinates to the header of the Qobj
        circuits_in_qobj = []
        indices = []
        theta_coordinates = []
        phi_coordinates = []
        while len(circuits_in_qobj) < max_experiments and index < num_pts:
            # until the length of circuits_in_qobj reaches the maximum allowed number of circuits to run, keep adding circuits
            circuits_in_qobj.append(circuits_transpiled[index])
            indices.append(index)
            theta_coordinates.append(target_points[index][0])
            phi_coordinates.append(target_points[index][1])
            index = index + 1
            # the loop has to break also if it reaches the index when a readout calibration has to be performed
            if index % readout_calibration_step == 0:
                break

        qobj_queue.put(assemble(circuits_in_qobj, backend=backend, shots=max_shots,
                                qobj_header=QobjHeader(is_readout=False, circuit_index=indices, theta=theta_coordinates,
                                                       phi=phi_coordinates)))

    # Run all the qobjs in the queue
    results_probabilities = []
    results_probabilities_corrected = []
    while not qobj_queue.empty():
        qobj = qobj_queue.get()
        # READOUT CALIBRATION
        if qobj.header.is_readout:
            readout_jobs = []
            print("Readout correction " + backend_identifier)
            if backend.provider() == QI:
                # If we get a readout circuit using QI, we know that we have to perform a readout calibration with 4
                # circuits
                for i in range(4):
                    print(i + 1, "/ 4")
                    readout_jobs.append(backend.run(qobj))  # Save the run job
                    readout_jobs[-1].result()  # Wait for the job to finish
                    qobj = qobj_queue.get()  # Get the next job from the queue
                # When the readout circuits are done, run the following job
            else:
                # If we are using IBM, there is only one readout job
                readout_jobs = backend.run(qobj)
                qobj = qobj_queue.get()  # Get the next job from the queue

            calculate_readout = True

        # Run the next qobj
        running_job = backend.run(qobj)
        if not backend.provider() == QI:
            # For IBM, wait only now for the results of the readout.
            # Doing so, we always send the readout circuits paired with the following batch of circuits
            readout_jobs.result()
            readout_jobs = [readout_jobs]

        if calculate_readout:  # Calculate readout only if necessary, not at every loop
            # Now it is possible to read the results and get the readout correction parameters
            readout_params = readout_analyze_data(readout_jobs, max_shots)  # Calculate readout parameters
            save_readout_parameters(path, backend_identifier, readout_params, index_copy1_transpiled,
                                    index_copy2_transpiled)
            calculate_readout = False

        # Wait for the results of the batch of cloning circuits
        running_job.result()
        # Calculate the marginal probabilities for the experiment that has just finished running
        results_probabilities_batch = analyze_data(running_job, index_copy1, index_copy2, max_shots)

        for i in range(len(results_probabilities_batch)):
            # Correct the results
            result_probabilities_corrected = correct_probabilities(results_probabilities_batch[i], readout_params)
            # Save both the results in a file
            save_result(qobj.header.theta[i], qobj.header.phi[i], results_probabilities_batch[i],
                        result_probabilities_corrected,
                        backend_identifier, path)
            # Append the single result to a list of all the results
            results_probabilities.append(results_probabilities_batch[i])
            results_probabilities_corrected.append(result_probabilities_corrected)

            # Print results
            print("========= " + backend_identifier + " ============")
            print(len(results_probabilities), "/", num_pts)
            print("Coordinates theta,phi: ", [qobj.header.theta[i], qobj.header.phi[i]])
            print("Fidelities copy1,copy2: ", [results_probabilities_batch[i][0], results_probabilities_batch[i][2]])
            print("Average fidelity copy1, copy2: ", [np.average(np.array(results_probabilities)[:, 0]),
                                                      np.average(np.array(results_probabilities)[:, 2])])
            print("========= corrected results ============")
            print(len(results_probabilities_corrected), "/", num_pts)
            print("Coordinates theta,phi: ", [qobj.header.theta[i], qobj.header.phi[i]])
            print("Fidelities copy1,copy2: ", [result_probabilities_corrected[0], result_probabilities_corrected[2]])
            print("Average fidelity copy1, copy2: ", [np.average(np.array(results_probabilities_corrected)[:, 0]),
                                                      np.average(np.array(results_probabilities_corrected)[:, 2])])
            print("=============================================")

    # Print the final results

    print("========= " + backend_identifier + " ============")
    print("FINISHED")
    print("===========================================")
    print("Final average fidelity copy1, copy2: ", [np.average(np.array(results_probabilities)[:, 0]),
                                                    np.average(np.array(results_probabilities)[:, 2])])
    print("===========================================")
    print("Final average fidelity copy1, copy2 with readout correction: ",
          [np.average(np.array(results_probabilities_corrected)[:, 0]),
           np.average(np.array(results_probabilities_corrected)[:, 2])])

    # Save final results on a file
    write_average_fidelity(path, backend_identifier, results_probabilities)
    write_average_fidelity(path, backend_identifier, results_probabilities_corrected, 'corrected')


if __name__ == '__main__':
    main()
