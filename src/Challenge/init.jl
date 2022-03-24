using QuantumCircuits
using QuantumCircuits.QML
using QuantumCircuits.QML.Optimization
using QuantumCircuits.QCircuits.Qiskit
using QuantumCircuits.QCircuits.Qiskit: qiskit
using QuantumCircuits.QCircuits.Circuit
using QuantumCircuits.QCircuits.Circuit: toQiskit, getCode
using QuantumCircuits.QCircuits.Math
using QuantumCircuits.Execute
using QuantumCircuits.Execute: setAndConvert, state2probability
using OpenScience21
using OpenScience21.Simulation.Gates

################################################################################

#using PyCall
qiskit.IBMQ.load_account()
provider = qiskit.IBMQ.get_provider(hub="ibm-q-community", group="ibmquantumawards", project="open-science-22")

jakarta = provider.get_backend("ibmq_jakarta")
# Simulated backend based on ibmq_jakarta's device noise profile
sim_noisy_jakarta = qiskit.providers.aer.QasmSimulator.from_backend(provider.get_backend("ibmq_jakarta"))
# Noiseless simulated backend
sim = qiskit.providers.aer.QasmSimulator()


################################################################################

const backend = QuantumSimulator()
const qiskitBackendSim = QiskitQuantum()
const qiskitBackendNoisySim = QiskitQuantum(sim_noisy_jakarta)
const qiskitBackendJakarta = QiskitQuantum(jakarta)

#err(x, y) = sum((abs.(x) .^ 2 - abs.(y) .^ 2) .^ 2) #/ length(x)
err(x, y) = sum(abs.(x - y))
err_state(x, y) = err(abs.(x) .^ 2, abs.(y) .^ 2)

################################################################################


function generate_circuit(trotter_steps, run_step, t=π, params=nothing, params2=nothing, params3=nothing; init=false, debug=false)
    qr = QuantumRegister(7, "q")
    qc = QCircuit(qr)

    # Prepare initial state (remember we are only evolving 3 of the 7 qubits on jakarta qubits (q_5, q_3, q_1) corresponding to the state |110>)
    if init
        qc.x([3, 5])  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
    end

    if params == nothing
        params = findU4paramsZZYYXX(t / trotter_steps, debug=debug)
        params2 = findU4paramsZZYYXX(t / (2 * trotter_steps), debug=debug)
        params3 = findU4paramsZZYYXXx2(t / (2 * trotter_steps), debug=debug)
    end

    qubits = [qr[1], qr[3], qr[5]]
    for s in 1:run_step
        isFirst = s == 1
        isLast = s == run_step

        trotter2U4(qc, qubits, t / trotter_steps, isFirst, isLast, params, params2, params3)
    end

    qc = decompose(qc)
    return qc, params, params2, params3
end

# function generate_circuit(trotter_steps, run_step, t=π, params=nothing; init=false, debug=false)
#     qr = QuantumRegister(7, "q")
#     qc = QCircuit(qr)
#
#     # Prepare initial state (remember we are only evolving 3 of the 7 qubits on jakarta qubits (q_5, q_3, q_1) corresponding to the state |110>)
#     if init
#         qc.x([3, 5])  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
#     end
#
#     if params == nothing
#         params = findU4paramsZZYYXX(t / trotter_steps, debug=debug)
#     end
#
#     qubits = [qr[1], qr[3], qr[5]]
#     for s in 1:run_step
#         trotterU4(qc, qubits, params)
#     end
#
#     qc = decompose(qc)
#     return qc, params
# end

function generate_empty_circuit(;init=false)
    qr = QuantumRegister(7, "q")
    qc = QCircuit(qr)

    # Prepare initial state (remember we are only evolving 3 of the 7 qubits on jakarta qubits (q_5, q_3, q_1) corresponding to the state |110>)
    if init
        qc.x([3, 5])  # DO NOT MODIFY (|q_5,q_3,q_1> = |110>)
    end

    return qc
end

################################################################################


using QuantumCircuits.QCircuits.Gates: Xmatrix, Ymatrix, Zmatrix

XXs = kron(kron(eye(2), Xmatrix), Xmatrix) + kron(kron(Xmatrix, Xmatrix), eye(2))
YYs = kron(kron(eye(2), Ymatrix), Ymatrix) + kron(kron(Ymatrix, Ymatrix), eye(2))
ZZs = kron(kron(eye(2), Zmatrix), Zmatrix) + kron(kron(Zmatrix, Zmatrix), eye(2))
Hs = XXs + YYs + ZZs

U_heis3(t) = exp(-im * Hs * t)

###################################################################################

function generate_ansact()
    n = 3
    qr = QuantumRegister(7, "q")
    qc = QCircuit(qr)
    qr = [qr[1], qr[3], qr[5]]

    qc.u3(qr)
    for i in (n-2):-1:0
        i = i+1
        qc.cx(qr[i], qr[i+1])
        #qc.rzx(qr[i], qr[i+1])
    end
    qc.u3(qr)
    for i in 0:(n-2)
        i = i+1
        qc.cx(qr[i], qr[i+1])
        #qc.rzx(qr[i], qr[i+1])
    end
    qc.u3(qr)

    return qc
end

################################################################################

function addMeasuresOS(qc)
    cr = ClassicalRegister(3)
    setClassicalRegister!(qc, cr)
    qc.measure([1, 3, 5], [0, 1, 2])
    nothing
end

################################################################################

function check_simulation_err(qc, t)
    sym_full = execute(backend, qc)
    exp_full = U_heis3(t) * ket"110"
    exp_full = abs.(exp_full) .^ 2
    err(sym_full, exp_full)
end

function create_full_circuit(qc_add, n)
    qr = QuantumRegister(7, "q")
    qc = QCircuit(qr)
    qc.x([3, 5])
    for i in 1:n
        append!(qc, qc_add)
    end
    addMeasuresOS(qc)

    return qc
end

################################################################################

function prepare_tomography_qc(qc, backend=nothing)
    qqc = toQiskit(qc)

    if backend != nothing
        oqc = optimize_pulses(qqc.qc, jakarta)
    else
        oqc = qqc.qc
    end

    qqr = getQRegister(qqc, "q")
    chec_qubits = [get(qqr, 1), get(qqr, 3), get(qqr, 5)]
    nqqc = qiskit.ignis.verification.tomography.state_tomography_circuits(oqc, chec_qubits)

    return nqqc
end

include("eval.jl")
