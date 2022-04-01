include("init.jl")

# state_end_expected = ket"110"
# state_start_sym = ket"110"
# state_start_sym_big = ket"0101000"
# state_zero = ket"000"

################################################################################


function generate_ansact_param(step_qc_exp, maxItr=1000)
    # Generate step qc
    step_qc = generate_empty_circuit()
    ansact = generate_ansact()
    # append
    append!(step_qc, step_qc_exp)
    append!(step_qc, ansact)
    length(getparameters(ansact))
    # Add measures
    addMeasuresOS(step_qc)

    # Set random parameters
    start_params = getRandParameters(step_qc)
    setparameters!(ansact, start_params)


    loss(params) = loss_expected_zero_state(execute(backend, step_qc, params))
    dloss(params) = real(loss'(params))
    of = OptimizationFunction(false, (x) -> (loss(x), dloss(x)), loss)

    val, xparams, itr = gradientDescent(of, start_params, α=0.01, maxItr=maxItr,
                                  argsArePeriodic=true, isExpectedZero=true, ϵ=1e-4, debug=false, useBigValInc=true)

    return xparams
end

trotter_steps = 10
trotter_steps_arr = [2, 2, 2, 2, 2]

# trotter_steps = 12
# trotter_steps_arr = [4, 2, 2, 2, 2]

# trotter_steps = 14
# trotter_steps_arr = [4, 2, 2, 2, 2, 2]

errs = []
errsai = []
ts = [π/4, π/2, 3*π/4, π]
#ts = [π]
best_param = []
best_qc = nothing
for t in ts
    println("Start for t=$t")
    qc, params, params2, params3 = generate_circuit(trotter_steps, trotter_steps, t)

    full_qc = generate_empty_circuit()
    full_qc.x([3, 5])
    append!(full_qc, qc)
    addMeasuresOS(full_qc)

    # caclulate the end state
    uerr = check_simulation_err(full_qc, t)

    ###########################################################################
    do_iter = 0
    do_check = true
    first_step = true
    inv_ansact = []
    for st in trotter_steps_arr
        step_qc_exp, _, _, _ = generate_circuit(trotter_steps, st, t, params, params2, params3, init=first_step)
        #step_qc_expmat = tomatrix(step_qc)
        do_iter += st
        println("Start AI TS $do_iter")

        if first_step
            opt_step_qc_exp = step_qc_exp
        else
            opt_step_qc_exp = generate_empty_circuit()
            append!(opt_step_qc_exp, inv_ansact)
            append!(opt_step_qc_exp, step_qc_exp)
        end

        # find best parameters
        best_params = generate_ansact_param(opt_step_qc_exp, 10000)

        # Generate ansact
        ansact = generate_ansact()
        setparameters!(ansact, best_params)
        inv_ansact = inv(ansact)
        bindparameters!(inv_ansact)

        # Check step reaults
        if do_check && do_iter < trotter_steps
            check_qc, _, _, _ = generate_circuit(trotter_steps, trotter_steps-do_iter, t, params, params2, params3)

            # Simulate full circuit
            check_step_qc_full = generate_empty_circuit()
            append!(check_step_qc_full, inv_ansact)
            append!(check_step_qc_full, check_qc)
            # Add measures
            addMeasuresOS(check_step_qc_full)

            # do check
            check_err = check_simulation_err(check_step_qc_full, t)
            check_err2 = check_circuits_err(full_qc, check_step_qc_full)
            println("Check value: $check_err, $check_err2")
            @assert check_err2 < 10e-2
        end

        # next initmat
        first_step = false
    end

    ###########################################################################
    addMeasuresOS(inv_ansact)
    check_err = check_simulation_err(inv_ansact, t)
    check_err2 = check_circuits_err(full_qc, inv_ansact)
    append!(errs, check_err)
    append!(errsai, check_err2)
end

println("==============================")
#println("The error is $(sum(errs)).")
println("The error Ai is $(mean(errsai)).")
println("==============================")
