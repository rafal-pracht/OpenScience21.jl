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
check_qc, _, _, _ = generate_circuit(trotter_steps, 2, t, params, params2, params3, init=false)
trotter_step_qc_exp, _, _, _ = generate_circuit(trotter_steps, st, t, params, params2, params3, init=false)
#### Step 1 - circuit ###############################
step3_params = [5.189896155542379, 4.707110397329062, 4.325905810068413, -0.8159924212917966, 3.110011220351916, 4.667709359716189, 1.0215152198560669, 5.116208100394577, 4.395063842373562, 4.370201901371649, 4.9135180721869425, 4.949975295863072, 4.637313944153916, 4.042830477313872, -0.6901894825683017, 1.1827639351923742, 3.0974447481049743, 4.022975117724089, 0.705887140029585, 3.514630564645274, 1.755049621148228, 4.151136862582433, 5.722357530716369, 0.8667811089902828, 0.3630963695624429, 0.1785816271535451, 3.4571312587872733]

tmp_ansact = generate_ansact()
step3_ansact = generate_ansact()
setparameters!(step3_ansact, step3_params)
step3_inv_ansact = inv(step3_ansact)
bindparameters!(step3_inv_ansact)
#### expected circiut ###############################
step_qc_exp = generate_empty_circuit()
append!(step_qc_exp, step3_inv_ansact)
append!(step_qc_exp, trotter_step_qc_exp)
#addMeasuresOS(step_qc_exp)

@assert length(getparameters(step_qc_exp)) == 0
err_state(tomatrix(check_qc) * tomatrix(step_qc_exp) * ket"0000000", ket"0101000")
err_state(tomatrix(check_qc) * tomatrix(step_qc_exp) * ket"0000000", tomatrix(qc_full) * ket"0000000")


tmp_qc_full = generate_empty_circuit()
append!(tmp_qc_full, step_qc_exp)
append!(tmp_qc_full, check_qc)
addMeasuresOS(tmp_qc_full)
execute(backend, tmp_qc_full)


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
# start_params = [1.5764568696923624, 1.0908808800022, 2.4252371170998916, 1.8639443060600367, 3.446628488779006, 4.377317661394165, 5.4941323482800914, 2.7412615189893237, 4.987825305684471, 3.2176798107893427, 1.9181246652054074, 6.103051613212941, 1.6943459153839973, 2.5283504304849904, 2.471571320001714, 0.363384700223248, 0.44494953380421215, 4.369514408776615, 5.8572397176591995, 4.9904251736390135, 2.6960480857237448, 3.7507514652393024, 2.9964093356989254, 2.8472290459219005, 2.0308048095101277, 0.5088982565685574, 2.288052006019751]
# start_params = [1.0808481661972744, 0.5669899130506986, 2.6451705381988964, 1.7404459205330207, 3.3052727253525407, 4.20764312227163, 4.937253903284112, 2.759119838259839, 4.793318765955065, 1.7409242381835606, 2.190273972584632, 5.501548807680859, 1.362408981041487, 2.2631548494399025, 2.0053070679102785, 0.12032411654173934, 0.4174475379244684, 4.644407333904697, 5.709216998663831, 5.097047763334831, 2.690316626109678, 3.9486279500908976, 2.9921660624519455, 2.5798040459155867, 1.9759784406454106, 0.5171690457217639, 2.3014311029811174]
# start_params = [1.1183704309752371, 0.5330049868725026, 2.5864863160919995, 1.517709635105656, 3.1379773116326626, 4.228131399915043, 4.839657602034448, 2.9578409509689414, 4.768586475194102, 1.6545055759000162, 2.9731956900142777, 5.4388456265418155, 1.3568252068474207, 2.2178334848189665, 2.016681314453245, -0.1393251937881101, 0.4175359667911205, 4.311315074278955, 5.654607248203027, 5.107649383590564, 2.6110369692045516, 3.964787457591229, 2.9897952619202135, 2.5850184346833758, 1.9933133206256328, 0.5174062767261149, 2.1454525723305404]
# start_params = [1.118701248839745, 0.5330109919755514, 2.5864455363659804, 1.5178030076960092, 3.1379811887238005, 4.228052136469158, 4.83951473643009, 2.9578555938528464, 4.768575419062368, 1.6545872882443375, 2.9732138548565903, 5.438798234834028, 1.3567873002796758, 2.2177969820984544, 2.0165985514741585, -0.1463770859679176, 0.4175308481100095, 4.311417477517902, 5.653628438550686, 5.107810333861132, 2.6110316764483468, 3.9648461032241777, 2.989792879025885, 2.585111455146589, 1.9933429040063286, 0.5174005280019287, 2.1454303587742265]
start_params = [1.1187171575015713, 0.5330092890891549, 2.586441510958004, 1.5178154678919513, 3.137963883844627, 4.228048049765462, 4.839509162353759, 2.9578752785288516, 4.768574523872119, 1.6545920069530968, 2.9732151650981935, 5.438796728232314, 1.3567895934397545, 2.2177922918734536, 2.0165705795863307, -0.14637716576453905, 0.4175300959136594, 4.311416597020989, 5.65363080884547, 5.1078093853571485, 2.6110312065705705, 3.9648509676298693, 2.989797720639175, 2.5851197335780856, 1.9933603087957705, 0.5174013121704684, 2.1454323563784383]

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

val, xparams, itr = gradientDescent(of, start_params, α=0.001, maxItr=20,
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
#append!(tmp_qc_full, step_qc_exp)
append!(tmp_qc_full, check_qc)
addMeasuresOS(tmp_qc_full)
execute(backend, tmp_qc_full)

##############################################################
start_params = [1.1187171575015713, 0.5330092890891549, 2.586441510958004, 1.5178154678919513, 3.137963883844627, 4.228048049765462, 4.839509162353759, 2.9578752785288516, 4.768574523872119, 1.6545920069530968, 2.9732151650981935, 5.438796728232314, 1.3567895934397545, 2.2177922918734536, 2.0165705795863307, -0.14637716576453905, 0.4175300959136594, 4.311416597020989, 5.65363080884547, 5.1078093853571485, 2.6110312065705705, 3.9648509676298693, 2.989797720639175, 2.5851197335780856, 1.9933603087957705, 0.5174013121704684, 2.1454323563784383]
tmp_ansact = generate_ansact()
setparameters!(tmp_ansact, start_params)

tmp_qc_full = generate_empty_circuit()
append!(tmp_qc_full, inv(tmp_ansact))
#append!(tmp_qc_full, step_qc_exp)
append!(tmp_qc_full, check_qc)
evaluate_qc_step4 = prepare_tomography_qc(tmp_qc_full)


jobs_step4 = evaluate_jobs(evaluate_qc_step4, sim_noisy_jakarta)
#jobs_step4 = evaluate_jobs(evaluate_qc_step4, jakarta)
# Job ID: 6220db43ff4b7c1bb06383b2
# Job ID: 6220db47d155c8af00e66721
# Job ID: 6220db4bd155c8f5c0e66722
# Job ID: 6220db4e5e625edf9babad54
# Job ID: 6220db523a4a26462ffbdbfc
# Job ID: 6220db563a4a2682c2fbdbfd
# Job ID: 6220db59d155c8bd3ce66723
# Job ID: 6220db5d3a4a262949fbdbfe
jobsid_step4 = ["6220db43ff4b7c1bb06383b2", "6220db47d155c8af00e66721", "6220db4bd155c8f5c0e66722",
                "6220db4e5e625edf9babad54", "6220db523a4a26462ffbdbfc", "6220db563a4a2682c2fbdbfd",
                "6220db59d155c8bd3ce66723", "6220db5d3a4a262949fbdbfe"]
jobs_step4 = [jakarta.retrieve_job(id) for id in jobsid_step4]


# Noisy sim: 0.639
# Real hardware: 0.490, 0.00479
results_step4 = evaluate_results(jobs_step4, evaluate_qc_step4)
