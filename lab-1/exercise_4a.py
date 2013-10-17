import matplotlib.pyplot as plt
import growth

def cobb_douglas_output(t, k, params):
    """
    Cobb-Douglas production function.

    Arguments:

        t:      (array-like) Time.
        k:      (array-like) Capital (per person/effective person).
        params: (dict) Dictionary of parameter values.
       
    Returns:

        y: (array-like) Output (per person/ effective person)

    """
    # extract params
    alpha = params['alpha']
    
    # Cobb-Douglas technology
    y = k**alpha
    
    return y
    
def marginal_product_capital(t, k, params):
    """
    Marginal product of capital with Cobb-Douglas production technology.

    Arguments:

        t:      (array-like) Time.
        k:      (array-like) Capital (per person/effective person).
        params: (dict) Dictionary of parameter values.
       
    Returns:

        y_k: (array-like) Derivative of output with respect to capital, k.

    """
    # extract params
    alpha = params['alpha']

    return alpha * k**(alpha - 1)

def equation_of_motion_capital(t, k, params):
    """
    Equation of motion for capital (per worker/effective worker).

    Arguments:

        t:      (array-like) Time.
        k:      (array-like) Capital (per person/effective person).
        params: (dict) Dictionary of parameter values.
       
    Returns:

        k_dot: (array-like) Time-derivative of capital (per worker/effective 
               worker).

    """
    # extract params
    s     = params['s']
    n     = params['n']
    g     = params['g']
    delta = params['delta']

    k_dot = s * cobb_douglas_output(t, k, params) - (n + g + delta) * k
    
    return k_dot
    
def solow_jacobian(t, k, params):
    """
    The Jacobian of the Solow model.
    
    Arguments:

        t:      (array-like) Time.
        k:      (array-like) Capital (per person/effective person).
        params: (dict) Dictionary of parameter values.
       
    Returns:

        jac: (array-like) Value of the derivative of the equation of 
             motion for capital (per worker/effective worker) with 
             respect to k.

    """
    # extract params
    s     = params['s']
    n     = params['n']
    g     = params['g']
    delta = params['delta']

    k_dot = s * marginal_product_capital(t, k, params) - (n + g + delta)
    
    return k_dot
    
def analytic_k_star(params): 
    """
    The steady-state level of capital stock per effective worker, k_bar, 
    in the Solow model is a function of the 5 exogenous parameters!
    
    """
    # extract params
    s     = params['s']
    n     = params['n']
    g     = params['g']
    alpha = params['alpha']
    delta = params['delta']
    
    return (s / (n + g + delta))**(1 / (1 - alpha))
     
# create a new model object
model = growth.SolowModel(cobb_douglas_output, marginal_product_capital, 
                          equation_of_motion_capital, solow_jacobian)

# create a dictionary of steady state expressions
steady_state_funcs = {'k_star':analytic_k_star}

# pass it as an argument to the set_steady_state_functions method
model.steady_state.set_functions(steady_state_funcs)

# calibrate the model and compute steady state values
growth.calibrate_cobb_douglas(model, 'GBR')

# create a new figure
fig = plt.figure(figsize=(12,8))

# irf for shock to alpha
model.plot_impulse_response('alpha', 1.5, 100, 'efficiency_units', False, True)

# display the figure
plt.show()

