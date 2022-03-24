include("../init.jl")


state_end_expected = ket"110"
state_start_sym = ket"110"

t = π
trotter_steps = 10
qc_full, params, params2, params3 = generate_circuit(trotter_steps, trotter_steps, t, init=true)
addMeasuresOS(qc_full)
# Execute simulation
check_simulation_err(qc_full, t)
#execute(backend, qc_full)


###############################

# Iterative algoritms
st = 2
check_qc, _, _, _ = generate_circuit(trotter_steps, 8, t, params, params2, params3, init=false)
step_qc_exp, _, _, _ = generate_circuit(trotter_steps, st, t, params, params2, params3, init=true)
@assert length(getparameters(step_qc_exp)) == 0
err(tomatrix(check_qc) * tomatrix(step_qc_exp) * ket"0000000", ket"0101000")
err(tomatrix(check_qc) * tomatrix(step_qc_exp) * ket"0000000", tomatrix(qc_full) * ket"0000000")

###############################

function getCheckValue(params)
    ansact = generate_ansact()
    setparameters!(ansact, params)
    inv_ansact = inv(ansact)

    check_step_qc_full = generate_empty_circuit()
    append!(check_step_qc_full, inv_ansact)
    append!(check_step_qc_full, check_qc)
    # Add measures
    addMeasuresOS(check_step_qc_full)
    check_simulation_err(check_step_qc_full, t)
end

function check_unitary_error(params)
    qc_tmp1 = generate_empty_circuit()
    append!(qc_tmp1, step_qc_exp)
    addMeasuresOS(qc_tmp1)
    tmp1_sv = qc_tmp1.measures_matrix * tomatrix(qc_tmp1) * ket"0000000"

    qc_tmp2 = generate_ansact()
    setparameters!(qc_tmp2, params)
    qc_tmp2 = inv(qc_tmp2)
    tmp2_sv = qc_tmp1.measures_matrix * tomatrix(qc_tmp2) * ket"0000000"

    return min_observe_unitary_error(tmp1_sv, tmp2_sv)
end


###############################

# Parameters
#start_params = [6.436585349626435, 3.5064762912929366, 5.417378549523721, 4.951710782201103, 2.8202725212610664, 2.6610028152685006, 1.493151979225325, 4.620028317366616, 3.3375281632180913, 1.942850689166652, -1.1747287791623802, 1.0819027892948092, 1.8327821316039612, 1.2347672501355669, 5.915790331657532, 5.457236741728066, 0.32787526669377787, -6.845344357522294, 4.857550796182953, 3.0147620689736554, 4.478928851291168, 1.3117161643159705, 1.9727355550126926, 3.6770740584075003, 5.385675638066927, 1.862252318664312, 5.053118574996168]
start_params = [6.490463894513841, 3.434943430067894, 5.375354561541277, 4.7906997798464115, 2.849325024133124, 2.7673620360366824, 1.4940499436354768, 4.487752591524359, 3.2948988322626382, 2.0225225509500486, -1.2343170007923738, 1.0431091490261049, 1.8837296276360387, 1.2169677999557502, 5.940807695305212, 5.421210358780618, 0.3173688878579456, -6.965567540712155, 4.835967484868803, 3.0127993733543517, 4.4620544822135475, 1.378337396849606, 1.980959799643362, 3.6789144092845376, 5.417860886313657, 1.8687391611628605, 5.092938820960128]
#start_params = getRandParameters(ansact)
getCheckValue(start_params)
check_unitary_error(start_params)


# Generate step qc
step_qc = generate_empty_circuit()
ansact = generate_ansact()
setparameters!(ansact, start_params)
# append
append!(step_qc, step_qc_exp)
append!(step_qc, ansact)
length(getparameters(ansact))
# Add measures
addMeasuresOS(step_qc)


#####################
loss(params) = loss_expected_zero_state(execute(backend, step_qc, params))
dloss(params) = real(loss'(params))
of = OptimizationFunction(false, (x) -> (loss(x), dloss(x)), loss)

of = OptimizationFunction(
        false,
        (x) -> qderivative(qiskitBackendNoisySim, step_qc, loss_expected_zero_state, x),
        (x) -> loss_expected_zero_state(execute(qiskitBackendNoisySim, setAndConvert(step_qc, x))))

of = OptimizationFunction(
        false,
        (x) -> qderivative(qiskitBackendJakarta, step_qc, loss_expected_zero_state, x),
        (x) -> loss_expected_zero_state(execute(qiskitBackendJakarta, setAndConvert(step_qc, x))))

#cp_start_params = start_params[:]
#start_params = cp_start_params[:]
#start_params = getRandParameters(ansact)
val, xparams, itr = gradientDescent(of, start_params, α=0.01, maxItr=50,
                              argsArePeriodic=true, isExpectedZero=true, ϵ=1e-8, debug=true, useBigValInc=true,
                              checkFn=(p) -> (getCheckValue(p), loss(p), check_unitary_error(p)))

#sum(abs.(xparams - start_params_new))
setparameters!(ansact, xparams)

##
setparameters!(step_qc, xparams)
start_params = xparams[:]
# start_params_74 = start_params_new[:]
# start_params_39 = start_params_new[:]
# start_params_36 = start_params_new[:]

println("The error is $val, $(loss(xparams))")
@assert abs(val - loss(xparams)) <= 1e-4 "Unexpected error in generate_ansact_param"
#return params, qc
#####################
#setparameters!(ansact, start_params)
inv_ansact = inv(ansact)


check_step_qc_full = generate_empty_circuit()
append!(check_step_qc_full, inv_ansact)
append!(check_step_qc_full, check_qc)
# Add measures
addMeasuresOS(check_step_qc_full)
check_simulation_err(check_step_qc_full, t)
execute(backend, check_step_qc_full)


################################################################################



err(tomatrix(check_step_qc_full) * ket"0000000", tomatrix(qc_full) * ket"0000000")

tomatrix(check_step_qc_full) * ket"0000000"



err(tomatrix(check_qc) * tomatrix(inv_ansact) * ket"0000000", ket"0101000")
err(qc_full.measures_matrix * tomatrix(check_qc) * tomatrix(inv_ansact) * ket"0000000", ket"110")
abs.(qc_full.measures_matrix * tomatrix(check_qc) * tomatrix(inv_ansact) * ket"0000000") .^ 2
abs.(qc_full.measures_matrix * tomatrix(qc_full) * ket"0000000") .^ 2

y, dy = qderivative(qiskitBackendNoisySim, step_qc, loss_expected_zero_state, xparams)
y2, dy2 = qderivative(qiskitBackendSim, step_qc, loss_expected_zero_state, xparams)

sum(abs.(dy - dy2))

check_ansact = generate_ansact()
setparameters!(check_ansact, xparams)
check_ansact = inv(check_ansact)
append!(check_ansact, check_qc)
cr = ClassicalRegister(3)
setClassicalRegister!(check_ansact, cr)
check_ansact.measure([1, 3, 5], [0, 1, 2])
check_full = execute(backend, qc_full)
err(sym_full, check_full)


#check_state = tomatrix(qc_full) * ket"0000000"
#check_state = tomatrix(check_qc) * tomatrix(step_qc) * ket"0000000"
#check_state = tomatrix(check_qc) * tomatrix(inv_ansact) * ket"0000000"



join(xparams, ", ")
"6.45462696829119, 3.4628539777662843, 5.381457724067881, 4.9087688572817365, 2.8362727497498863, 2.815284242011393, 1.4559848200439092, 4.557887265162635, 3.2191782562021074, 2.017749990259502, -1.2773778371627664, 1.038280475768157, 1.7763927041852872, 1.2117391094655119, 5.95751465521469, 5.455912091958163, 0.32787526669377787, -6.947993415522681, 4.87336831611616, 3.071494811643587, 4.542395308051386, 1.3814250769083867, 1.9727355550126926, 3.6810723669714087, 5.306243497567738, 1.862252318664312, 5.133659675982915"

of.f(xparams)
xparams



abs.(qc_full.measures_matrix * check_state)
for (i, v) in enumerate(abs.(check_state))
    if v > 0.01
        @show i, v
    end
end

out = Dict("000"=>5294, "001"=>315, "010"=>503, "011"=>236, "100"=>525, "101"=>278, "110"=>782, "111"=>259)
shots = sum([c for (k, c) in out])
# 0.27030695561468754
using QGen.Quantum.Execute: extractProbability
pro_mes = extractProbability(out, 3, shots)
loss_expected_zero_state(pro_mes)

execute(backend, step_qc, xparams)
execute(backend, step_qc, start_params)

loss_expected_zero_state(execute(backend, step_qc, xparams))
pro_mes = execute(qiskitBackendJakarta, setAndConvert(step_qc, xparams))

pro_mes2 = execute(qiskitBackendSim, setAndConvert(step_qc, xparams))
pro_mes3 = execute(qiskitBackendNoisySim, setAndConvert(step_qc, xparams))

execute(qiskitBackendSim, setAndConvert(step_qc, start_params))
execute(qiskitBackendNoisySim, setAndConvert(step_qc, start_params))

loss_expected_zero_state(pro_mes)


# Rand params
params = getRandParameters(ansact)
setparameters!(ansact, params)

setparameters!(step_qc, params)

49
start_param = [
3.10533603728926
3.112581816710326
2.6463145809437174
6.133225716695942
3.78636643176749
1.8187176425437817
0.375079482875157
0.6571266498955395
2.801217447829762
0.5916427818203986
0.5790649922699331
0.712208515948833
1.9966225293108957
1.9593341367998065
2.2290937876201835
1.1749078516564981
5.343155591107844
4.584391363834062
2.36219367045622
0.5612657644487052
1.1180274074966425
0.7111816437845817
5.861886540549156
2.4199006556390135
4.710833068188319
1.687246263047568
0.5850192294604584
]



per = rand(27) .- 0.5
param = start_param + per


4.046576639909423
2.1094424943285826
0.9201902642387239
4.770530249750385
4.7647761819852565
6.204241817619136
1.898538156438555
5.785115385976787
3.192476843851657
3.7461808924516475
5.18730667030047
5.107773209523165
3.2320049641622304
1.9056591555970601
3.9177329454116303
2.18314166981667
3.7858593406069185
2.1926428624140297
0.05231222345829235
6.190847811985151
5.868307355303737
5.944672797969161
0.8045669169611496
6.090000186546791
5.376862289442649
2.708270293885021
1.0663248795742288


st = 2
check_qc, _, _, _ = generate_circuit(trotter_steps, 8, t, params, params2, params3, init=false)
step_qc_exp, _, _, _ = generate_circuit(trotter_steps, st, t, params, params2, params3, init=true)
@assert length(getparameters(step_qc_exp)) == 0
err(tomatrix(check_qc) * tomatrix(step_qc_exp) * ket"0000000", ket"0101000")
err(tomatrix(check_qc) * tomatrix(step_qc_exp) * ket"0000000", tomatrix(qc_full) * ket"0000000")

###############################
