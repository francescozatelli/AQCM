from qiskit import IBMQ
from quantuminspire.api import QuantumInspireAPI
from quantuminspire.credentials import save_account
from quantuminspire.qiskit import QI

"""
Here we get the backend associated to a nickname that we decided for the backend for convenience
"""

TOKEN_QI = 'YOUR_TOKEN'
TOKEN_IBMQ = 'YOUR_TOKEN '


dict_backend_identifiers_ibm = {'simulator': 'ibmq_qasm_simulator',
                                'ibmqx2': 'ibmqx2',
                                'melbourne': 'ibmq_16_melbourne',
                                'vigo': 'ibmq_vigo',
                                'ourense': 'ibmq_ourense',
                                'valencia': 'ibmq_valencia',
                                'athens': 'ibmq_athens',
                                'santiago': 'ibmq_santiago'}

dict_backend_identifiers_qi = {'spin2': 'Spin-2',
                               'starmon5': 'Starmon-5',
                               'qi_simulator': 'QX Single-Node Simulator'}


def find_backend_by_nickname(backend_identifier):
    backend_nickname = str(backend_identifier)

    # Check if the backend nickname belongs to IBM and get it if that's the case
    if backend_nickname in dict_backend_identifiers_ibm:
        IBMQ.save_account(TOKEN_IBMQ, overwrite=True)
        provider = IBMQ.load_account()
        backend = provider.get_backend(dict_backend_identifiers_ibm[backend_nickname])
    # Check if the backend nickname belongs to QI and get it if that's the case
    elif backend_nickname in dict_backend_identifiers_qi:
        save_account(TOKEN_QI)
        QI.set_authentication()
        backend = QI.get_backend(dict_backend_identifiers_qi[backend_nickname])
    else:
        print('Invalid Backend!')
        exit()
    #Return the backend
    return backend