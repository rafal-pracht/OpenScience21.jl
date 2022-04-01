# Hi everyone!
# I hope your projects are going well, and you're having a fun time of it! :)
# I have a quick clarification for those planning to do some kind of
# zero-noise extrapolation. Be sure to not directly extrapolate the state
# tomography fidelity to zero noise. Extrapolating other expectation values
# that then inform data going into the state tomography calculation is just fine though.
#
# Stefan H. Sack  15 days ago
# If we mitigate the data that goes into the StateTomographyFitter,
# how can we get it to be the correct data type? The function only takes
# Result() objects and not just a counts dictionary.
#
# AJ Rasmusson  12 days ago
# There are a couple options, but they depend on the details of your solution.
# Let's chat in a DM
backend_sim_jakarta = QiskitQuantum(sim_noisy_jakarta)

using QuantumCircuits.QCircuits.Gates: CX

function generate_3cx_circuit(qc)
    newQC = generate_empty_circuit()
    for c in getCode(qc)
        if typeof(c) == CX
            add!(newQC, c)
            newQC.barrier()
            add!(newQC, c)
            newQC.barrier()
            add!(newQC, c)
        else
            add!(newQC, c)
        end
    end

    return newQC
end



qc = QCircuit(2)
qc.h(0)
qc.cx(0, 1)
qc.x(1)
qc.measure([0, 1], [0, 1])


qc2 = QCircuit(2)
qc2.h(0)
qc2.cx(0, 1)
qc2.barrier()
qc2.cx(0, 1)
qc2.barrier()
qc2.cx(0, 1)
qc2.barrier()
qc2.x(1)
qc2.measure([0, 1], [0, 1])

qc = tmp_qc_full
addMeasuresOS(qc)

M0 = execute(backend, qc)
M1 = execute(backend_sim_jakarta, toQiskit(qc))

qc2 = generate_3cx_circuit(qc)
addMeasuresOS(qc2)
M3 = execute(backend_sim_jakarta, toQiskit(qc2))

execute(backend, qc2)

calc_err(m1, m2) = sum(abs.(m1 - m2))
calc_err(M0, M1)
calc_err(M0, M3)
Mex = 3/2 * M1 - 1/2 * M3
calc_err(M0, Mex)
