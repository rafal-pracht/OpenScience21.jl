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
check_qc, _, _, _ = generate_circuit(trotter_steps, 4, t, params, params2, params3, init=false)
trotter_step_qc_exp, _, _, _ = generate_circuit(trotter_steps, st, t, params, params2, params3, init=false)
#### Step 1 - circuit ###############################
step2_params = [6.357739083193767, 4.110466988994183, 4.566133708100635, -4.727581603064405, 6.207391247986331, -1.3477711816915519, -0.7698439582879397, 4.677177128084215, -1.88364164245155, -6.961405873009646, 2.3088656324648755, 1.9854073372986245, -2.4954433452385736, 6.326779366664302, -0.7789648567030898, 2.995616831369772, 6.237335880828537, 0.6943207146536418, 0.37027914593259564, 0.9373275950100715, 2.0325673904575456, 1.3926281522358193, 5.042496816647515, 0.3933698184667445, 3.4025273194720507, 0.9766806621264883, 3.8089690519454855]

tmp_ansact = generate_ansact()
step2_ansact = generate_ansact()
setparameters!(step2_ansact, step2_params)
step2_inv_ansact = inv(step2_ansact)
bindparameters!(step2_inv_ansact)
#### expected circiut ###############################
step_qc_exp = generate_empty_circuit()
append!(step_qc_exp, step2_inv_ansact)
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
# start_params = [5.510624874968758, 4.731878064388753, 4.900181698295356, -0.7772941524378829, 4.288753243775524, 4.534026947674572, 0.9628921021029925, 5.052499307052674, 4.3163477736916755, 4.090039770917148, 5.146232943849546, 4.834904355385483, 4.298421076980501, 4.298947426088873, -0.2245520107841166, 1.8269035949274577, 3.095670931343558, 4.2964291187218375, 0.6880915325673521, 2.8018895475710592, 1.755438189888574, 3.6909076939142063, 5.721085576526464, 7.378200980981615, 0.41420998486446187, 0.17120739482060593, 10.088782846749732]
# start_params = [5.187965455615748, 4.7072210399028345, 4.326505144053168, -0.8158941350039923, 3.111980277303591, 4.665345765167737, 1.0270995070313362, 5.088159589626752, 4.398365594174959, 4.3705576878782875, 4.917549761672741, 4.950009003101304, 4.629179967940117, 4.050156889347045, -0.6840998229626889, 1.1977756653867702, 3.097444756172329, 4.023038389469699, 0.7060636965564819, 3.513475977641944, 1.7473365924622046, 4.156426317899939, 5.722349725139118, 7.139067891637897, 0.3706195248147545, 0.1785748671934303, 3.45711420264948]
# start_params = [5.189909571913178, 4.7061718348030634, 4.325920298601799, -0.8159940610859193, 3.1100111597039337, 4.667711202614713, 1.0215012685993325, 5.116278529322291, 4.39485224831755, 4.370199934115478, 4.912671371874754, 4.949975882952213, 4.63731372146988, 4.042827546675875, -0.689958025429904, 1.182609375802692, 3.0974449271310225, 4.022880155085733, 0.7058858557420989, 3.5146428736574116, 1.7550444429952092, 4.151140142069366, 5.722356721254833, 0.866787454432816, 0.3630921350881019, 0.1785746966481526, 3.4571302116970064]
start_params = [5.189896155542379, 4.707110397329062, 4.325905810068413, -0.8159924212917966, 3.110011220351916, 4.667709359716189, 1.0215152198560669, 5.116208100394577, 4.395063842373562, 4.370201901371649, 4.9135180721869425, 4.949975295863072, 4.637313944153916, 4.042830477313872, -0.6901894825683017, 1.1827639351923742, 3.0974447481049743, 4.022975117724089, 0.705887140029585, 3.514630564645274, 1.755049621148228, 4.151136862582433, 5.722357530716369, 0.8667811089902828, 0.3630963695624429, 0.1785816271535451, 3.4571312587872733]
getCheckValue(start_params)

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

#start_params = getRandParameters(ansact)

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
start_params = [5.189896155542379, 4.707110397329062, 4.325905810068413, -0.8159924212917966, 3.110011220351916, 4.667709359716189, 1.0215152198560669, 5.116208100394577, 4.395063842373562, 4.370201901371649, 4.9135180721869425, 4.949975295863072, 4.637313944153916, 4.042830477313872, -0.6901894825683017, 1.1827639351923742, 3.0974447481049743, 4.022975117724089, 0.705887140029585, 3.514630564645274, 1.755049621148228, 4.151136862582433, 5.722357530716369, 0.8667811089902828, 0.3630963695624429, 0.1785816271535451, 3.4571312587872733]
tmp_ansact = generate_ansact()
setparameters!(tmp_ansact, start_params)

tmp_qc_full = generate_empty_circuit()
append!(tmp_qc_full, inv(tmp_ansact))
#append!(tmp_qc_full, step_qc_exp)
append!(tmp_qc_full, check_qc)
evaluate_qc_step3 = prepare_tomography_qc(tmp_qc_full)


jobs_step3 = evaluate_jobs(evaluate_qc_step3, sim_noisy_jakarta)
#jobs_step3 = evaluate_jobs(evaluate_qc_step3, jakarta)
# Job ID: 621fa4d4e71dc7981752fb33
# Job ID: 621fa4d9e71dc7920b52fb34
# Job ID: 621fa4de5e625e8154abaa4f
# Job ID: 621fa4e35e625ed478abaa50
# Job ID: 621fa4e64668ab31598d9efd
# Job ID: 621fa4eb4668ab2b898d9efe
# Job ID: 621fa4ef3e2bc416aa4bf979
# Job ID: 621fa4f33a4a262632fbd90c
jobsid_step3 = ["621fa4d4e71dc7981752fb33", "621fa4d9e71dc7920b52fb34", "621fa4de5e625e8154abaa4f",
                "621fa4e35e625ed478abaa50", "621fa4e64668ab31598d9efd", "621fa4eb4668ab2b898d9efe",
                "621fa4ef3e2bc416aa4bf979", "621fa4f33a4a262632fbd90c"]
jobs_step3 = [jakarta.retrieve_job(id) for id in jobsid_step3]

# Noisy sim: 0.517
# Real hardware: 0.38, 0.173
results_step3 = evaluate_results(jobs_step3, evaluate_qc_step3)
