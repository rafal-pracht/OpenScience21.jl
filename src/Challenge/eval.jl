using PyCall
using Statistics

py"""
import numpy as np
import matplotlib.pyplot as plt

# Importing standard Qiskit modules
from qiskit import QuantumCircuit, QuantumRegister, IBMQ, execute, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.tools.monitor import job_monitor
from qiskit.circuit import Parameter

# Import state tomography modules
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
from qiskit.quantum_info import state_fidelity

# Import qubit states Zero (|0>) and One (|1>), and Pauli operators (X, Y, Z)
from qiskit.opflow import Zero, One, I, X, Y, Z


# Compute the state tomography based on the st_qcs quantum circuits and the results from those ciricuits
def state_tomo(result, st_qcs):
    # The expected final state; necessary to determine state tomography fidelity
    target_state = (One ^ One ^ Zero).to_matrix()  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
    # Fit state tomography results
    tomo_fitter = StateTomographyFitter(result, st_qcs)
    rho_fit = tomo_fitter.fit(method='lstsq')
    # Compute fidelity
    fid = state_fidelity(rho_fit, target_state)
    return fid


def evaluate(qc: QuantumCircuit, backend, shots=8192, reps=8):
    jobs = []
    for _ in range(reps):
        # execute
        job = execute(qc, backend, shots=shots)
        print('Job ID', job.job_id())
        jobs.append(job)

    # wait for finish
    for job in jobs:
        job_monitor(job)
        try:
            if job.error_message() is not None:
                print(job.error_message())
        except:
            pass

    # Compute tomography fidelities for each repetition
    fids = []
    for job in jobs:
        fid = state_tomo(job.result(), qc)
        fids.append(fid)

    return np.mean(fids), np.std(fids)
"""


function evaluate_jobs(qc, backend; shots=8192, reps=8)
    println("== evaluate_jobs ==")
    jobs = []
    for _ in 1:reps
        # execute
        job = qiskit.execute(qc, backend, shots=shots)
        println("Job ID: $(job.job_id())")
        push!(jobs, job)
    end

    return jobs
end

function evaluate_results(jobs, qc)
    # wait for finish
    for job in jobs
        qiskit.tools.monitor.job_monitor(job)
    end

    # Compute tomography fidelities for each repetition
    fids = []
    for job in jobs
        fid = py"state_tomo"(job.result(), qc)
        push!(fids, fid)
    end

    return mean(fids), std(fids)
end

function evaluate(qc, backend; shots=8192, reps=8)
    jobs = evaluate_jobs(qc, backend; shots=shots, reps=reps)

    return evaluate_results(jobs, qc)
end
