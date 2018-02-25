import numpy as np
import scipy.optimize

from optimization.poissonCalculator import ALMOST_ZERO


def optimize(config, C, A, initial_x, alphas, edge_log_probabilities_array):
    x = initial_x
    bounds = [(ALMOST_ZERO, 1)] * initial_x.shape[0]
    len_constraints = A.shape[0]
    edge_log_probabilities_array = np.array(edge_log_probabilities_array)

    lambd = np.zeros(len_constraints)
    k = 0
    mi = config.penalty_coefficient  # 0.00001
    last_x = np.nan
    last_goal_fun = np.nan
    last_error = np.nan
    results = []
    no_change_iterations = 0
    while k < config.max_iteration:
        print("")
        print("iteration: ", k)
        result = scipy.optimize.minimize(functionToOptimize, x, args=(C, A, edge_log_probabilities_array, alphas, mi, lambd), method='L-BFGS-B',
                                         jac=jacobian,
                                         bounds=bounds,
                                         # options = {"gtol" : 1e-7, "disp":True},
                                         options={"ftol": 1e-9,
                                                  # The iteration stops when (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1}
                                                  # <= ftol.
                                                  "gtol": 1e-9,
                                                  # The iteration will stop when max{|proj g_i | i = 1, ..., n}
                                                  # <= gtol where pg_i is the i-th component of the projected gradient.
                                                  "disp": False}
                                         )  # callback = callback_fun) #Doesn't work?!
        x = result.x
        goal_fun = goalFunction(x, edge_log_probabilities_array, alphas)
        error = sumOfQuadraticConstraints(constraints(x, C, A))
        goal_change = abs(goal_fun - last_goal_fun) / abs(last_goal_fun)
        error_change = abs(last_error - error) / last_error
        x_diff = np.sum(np.power(x - last_x, 2))
        print("result.fun={}, goalFunction={}, error={}".format(result.fun, goal_fun, error))
        if k > 0:
            print("goalFunction.change={}, xdiff.norm={}, error.change={}".format(goal_change,
                                                                                  x_diff, error_change))
        last_goal_fun = goal_fun
        last_x = x
        last_error = error
        results.append((goal_fun, error, goal_change, error_change, x_diff))
        if x_diff < config.min_x_change:
            no_change_iterations += 1
            if no_change_iterations == config.min_x_change_iterations:
                break
        else:
            no_change_iterations = 0
        print("----------")
        lambd = updateLambda(lambd, mi, x, C, A)
        mi = config.penalty_coefficient_multiplier * mi
        k += 1

    print("Optymization finished")

    return x, results


def constraints(x, C, A):
    return C@x-A


# goal function
def goalFunction(x, edge_log_probabilities_array, alphas):
    return edge_log_probabilities_array @ x + alphas @ np.log(x)


def der_goalFunction(x, edge_log_probabilities_array, alphas):
    return np.add(edge_log_probabilities_array, alphas / x)


def sumOfQuadraticConstraints(c):
    return np.sum(np.power(c, 2))


def sumOfConstraints(c, lambd):
    return c @ lambd


def updateLambda(lambd, mi, x, C, A):
    return lambd - constraints(x, C, A) * mi


def functionToOptimize(x, C, A, edge_log_probabilities_array, alphas, mi, lambd):
    c = constraints(x, C, A)
    # res =
    minGF = - goalFunction(x, edge_log_probabilities_array, alphas)
    miSumQC = + mi / 2 * sumOfQuadraticConstraints(c)
    minSC = - sumOfConstraints(c, lambd)
    res = minGF + miSumQC + minSC
    return res


def jacobian(x, C, A, edge_log_probabilities_array, alphas, mi, lambd):
    # return np.ones(x.shape)
    c = constraints(x, C, A)
    jac_goalFunction = - der_goalFunction(x, edge_log_probabilities_array, alphas)
    jac_quadraticConstraints = mi * c.T @ C
    jac_sumOfConstraints = - lambd @ C
    # if __DEBUG:
    #    print("x.dim={}".format(x.shape))
    #    print("jac_gF.dim={}".format(jac_goalFunction.shape))
    #    print("jac_qC.dim={}".format(jac_quadraticConstraints.shape))
    #    print("jac_soC.dim={}".format(jac_sumOfConstraints.shape))
    result = np.add(jac_goalFunction, np.add(jac_quadraticConstraints, jac_sumOfConstraints))
    return result