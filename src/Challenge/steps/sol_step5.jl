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
trotter_step_qc_exp, _, _, _ = generate_circuit(trotter_steps, st, t, params, params2, params3, init=false)
#### Step 1 - circuit ###############################
step4_params = [1.1187171575015713, 0.5330092890891549, 2.586441510958004, 1.5178154678919513, 3.137963883844627, 4.228048049765462, 4.839509162353759, 2.9578752785288516, 4.768574523872119, 1.6545920069530968, 2.9732151650981935, 5.438796728232314, 1.3567895934397545, 2.2177922918734536, 2.0165705795863307, -0.14637716576453905, 0.4175300959136594, 4.311416597020989, 5.65363080884547, 5.1078093853571485, 2.6110312065705705, 3.9648509676298693, 2.989797720639175, 2.5851197335780856, 1.9933603087957705, 0.5174013121704684, 2.1454323563784383]

tmp_ansact = generate_ansact()
step4_ansact = generate_ansact()
setparameters!(step4_ansact, step4_params)
step4_inv_ansact = inv(step4_ansact)
bindparameters!(step4_inv_ansact)
#### expected circiut ###############################
step_qc_exp = generate_empty_circuit()
append!(step_qc_exp, step4_inv_ansact)
append!(step_qc_exp, trotter_step_qc_exp)
#addMeasuresOS(step_qc_exp)

@assert length(getparameters(step_qc_exp)) == 0
err_state(tomatrix(step_qc_exp) * ket"0000000", ket"0101000")
err_state(tomatrix(step_qc_exp) * ket"0000000", tomatrix(qc_full) * ket"0000000")


tmp_qc_full = generate_empty_circuit()
append!(tmp_qc_full, step_qc_exp)
addMeasuresOS(tmp_qc_full)
execute(backend, tmp_qc_full)

###############################

function getCheckValue(params)
    ansact = generate_ansact()
    setparameters!(ansact, params)
    inv_ansact = inv(ansact)

    check_step_qc_full = generate_empty_circuit()
    append!(check_step_qc_full, inv_ansact)
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
# start_params = [1.3074139897354728, 0.5336490045385581, 2.4568535729154437, 3.318798857322259, 0.23154315417368757, 4.109186956787016, 3.1327696958680473, 6.196015337198003, 5.483788640849114, 0.07228688236343765, 6.063108055980849, 5.897844825526794, 4.99787283562687, 3.0203366040307267, 4.313080680970031, 0.1926720605196993, 2.2213780811581074, 1.7450434426558326, 4.940088416493364, 4.926796649597691, 5.2158253353061665, 5.9647454624186365, 1.4978941361061908, 6.136918858144318, 3.345066280621706, 1.808332924079646, 1.914145236424726]
start_params = [1.103633523141456, 0.5621723387373332, 3.535362761806341, 3.1465964893594798, 0.2416797028860025, 3.7901137174699655, 2.515512466095495, 6.762078821484598, 5.187907546418758, -0.19860568802167688, 5.997524383569426, 5.905306784345511, 6.0845429036054846, 3.1437960954072843, 4.530029341809641, 1.2337272277202151, 2.221508726842684, 1.6606762213540665, 4.12874467494519, 4.827888262840575, 5.603301739361324, 5.894694449973559, 1.4978959576376278, 6.261167575543967, 2.7695944642836734, 1.8088920099024433, 1.9134833740486716]

# Generate step qc
step_qc = generate_empty_circuit()
ansact = generate_ansact()
# append
append!(step_qc, step_qc_exp)
append!(step_qc, ansact)
@assert length(getparameters(ansact)) == 27
setparameters!(ansact, start_params)
# Add measures
addMeasuresOS(step_qc)
# start_params = getRandParameters(ansact)
# println(join(start_params, ", "))

#####################
loss(params) = loss_expected_zero_state(execute(backend, step_qc, params))
dloss(params) = real(loss'(params))
#of = OptimizationFunction(false, (x) -> (loss(x), dloss(x)), loss)

function call_N_qderivative(x, N)
    y  = 0.0
    dy = zeros(27)
    for i in 1:N
        tmp_y, tmp_dy = qderivative(qiskitBackendNoisySim, step_qc, loss_expected_zero_state, x)
        #tmp_y, tmp_dy = qderivative(qiskitBackendSim, step_qc, loss_expected_zero_state, x)
        y = y + tmp_y
        #println(tmp_y)
        dy = dy + tmp_dy
    end
    return y/N, dy/N
end

function call_N_loss_expected_zero_state(x, N)
    y  = 0.0
    for i in 1:N
        #tmp_y = loss_expected_zero_state(execute(qiskitBackendNoisySim, setAndConvert(step_qc, x)))
        tmp_y = qexecute(qiskitBackendNoisySim, step_qc, loss_expected_zero_state, x)
        y = y + tmp_y
        #println(tmp_y)
    end
    return y/N
end

of = OptimizationFunction(
        false,
        (x) -> call_N_qderivative(x, 30),
        (x) -> call_N_loss_expected_zero_state(x, 30))

##############################
getCheckValue(start_params)
loss(start_params)

val, xparams, itr = gradientDescent(of, start_params, α=0.1, maxItr=20,
                              argsArePeriodic=true, isExpectedZero=true, ϵ=1e-8, debug=true, useBigValInc=true,
                              checkFn=(p) -> (getCheckValue(p), loss(p), check_unitary_error(p)))

getCheckValue(xparams)
loss(xparams)

######### Check ####
tmp_ansact = generate_ansact()
#setparameters!(tmp_ansact, xparams)
setparameters!(tmp_ansact, start_params)

tmp_qc_full = generate_empty_circuit()
append!(tmp_qc_full, inv(tmp_ansact))
addMeasuresOS(tmp_qc_full)
execute(backend, tmp_qc_full)

##############################################################
start_params = [1.103633523141456, 0.5621723387373332, 3.535362761806341, 3.1465964893594798, 0.2416797028860025, 3.7901137174699655, 2.515512466095495, 6.762078821484598, 5.187907546418758, -0.19860568802167688, 5.997524383569426, 5.905306784345511, 6.0845429036054846, 3.1437960954072843, 4.530029341809641, 1.2337272277202151, 2.221508726842684, 1.6606762213540665, 4.12874467494519, 4.827888262840575, 5.603301739361324, 5.894694449973559, 1.4978959576376278, 6.261167575543967, 2.7695944642836734, 1.8088920099024433, 1.9134833740486716]
tmp_ansact = generate_ansact()
setparameters!(tmp_ansact, start_params)

tmp_qc_full = generate_empty_circuit()
append!(tmp_qc_full, inv(tmp_ansact))
evaluate_qc_step5 = prepare_tomography_qc(tmp_qc_full)

jobs_step5 = evaluate_jobs(evaluate_qc_step5, sim_noisy_jakarta)
#jobs_step5 = evaluate_jobs(evaluate_qc_step5, jakarta)
# Job ID: 6220f14f4668ab66ef8da277
# Job ID: 6220f153877263eede6fc339
# Job ID: 6220f157d155c810fde66770
# Job ID: 6220f1613a4a260a7ffbdc34
# Job ID: 6220f164e71dc766df52fe82
# Job ID: 6220f16fe71dc777cf52fe84
# Job ID: 6220f172e71dc7046f52fe85
# Job ID: 6220f174762eb758e181c8c6
jobsid_step5 = ["6220f14f4668ab66ef8da277", "6220f153877263eede6fc339", "6220f157d155c810fde66770",
                "6220f1613a4a260a7ffbdc34", "6220f164e71dc766df52fe82", "6220f16fe71dc777cf52fe84",
                "6220f172e71dc7046f52fe85", "6220f174762eb758e181c8c6"]

# Job ID: 6230d605d73e0a4b4a31ff7c
# Job ID: 6230d609d10f746d4467ae33
# Job ID: 6230d60c74a4dc7a8535b5ae
# Job ID: 6230d610d10f74540667ae34
# Job ID: 6230d6122001bb71d86e9827
# Job ID: 6230d6150d6e0d97d1138696
# Job ID: 6230d618314f79a95e20c20d
# Job ID: 6230d61b0d6e0d4cd6138697
jobsid_step5 = ["6230d605d73e0a4b4a31ff7c", "6230d609d10f746d4467ae33", "6230d60c74a4dc7a8535b5ae",
                "6230d610d10f74540667ae34", "6230d6122001bb71d86e9827", "6230d6150d6e0d97d1138696",
                "6230d618314f79a95e20c20d", "6230d61b0d6e0d4cd6138697"]


jobs_step5 = [jakarta.retrieve_job(id) for id in jobsid_step5]

# Noisy sim: 0.797, 0.826
# Real hardware:  0.783, 0.00448
# Real hardware2: 0.7975091601169371, 0.011794318866976355
results_step5 = evaluate_results(jobs_step5, evaluate_qc_step5)


# Slack
# Noisy sim   : ~0.832
# Real device : ~0.765
