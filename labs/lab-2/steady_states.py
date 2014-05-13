import sympy as sp
from scipy import linalg, optimize

class SteadyState(object):
    """Abstract class representing the deterministic steady state of a model."""

    def __init__(self, model):
        """
        Initializes a SteadyState object with the following attributes:
        
            model: (object) An instance of some Model class.
            
        """
        self.model        = model
        self.eigenvalues  = None
        self.eigenvectors = None
        
    def set_functions(self, func_dict):
        """
        Modifies the model's steady_state_functions attribute.
        
        Arguments:
            
            func_dict: (dict) Dictionary of analytic function defining the 
                       model's steady state as functions of model parameters.
                       
        """  
        self.functions = func_dict
        
    def set_values(self):
        """
        Computes the steady state values using the dictionaries of analytic
        functions and model parameters.
                  
        """
        # store steady state values in a dictionary
        steady_state_values = {}
        
        # populate the dictionary of steady state values
        for key, func in self.functions.iteritems():
            steady_state_values[key] = func(self.model.params)
        
        self.values = steady_state_values
    
    def find_values(self, func, method, **kwargs):
        """
        Provides an interface to the various methods for finding the root of 
        a univariate or multivariate function in scipy.optimize.
    
        Arguments:
            
            method:   (str) For univariate functions method must be one of: 
                      'brentq', 'brenth', 'ridder', 'bisect', or 'newton'; for
                      multivariate function method must be one of: 'hybr', 'lm',
                      'broyden1', 'broyden2', 'anderson', 'linearmixing', 
                      'diagbroyden', 'excitingmixing', 'krylov'.
            
            **kwargs: (dict) Dictionary of method specific keyword arguments.
        
        """ 
        # list of valid multivariate root-finding methods
        multivariate_methods = ['hybr', 'lm', 'broyden1', 'broyden2', 
                                'anderson', 'linearmixing', 'diagbroyden', 
                                'excitingmixing', 'krylov']  
        # univariate methods     
        if method == 'brentq':
            res = optimize.brentq(func, args=self.model.args, **kwargs)
        elif method == 'brenth':
            res = optimize.brenth(func, args=self.model.args, **kwargs)
        elif method == 'ridder':
            res = optimize.ridder(func, args=self.model.args, **kwargs)
        elif method == 'bisect':
            res = optimize.bisect(func, args=self.model.args, **kwargs)
        elif method == 'newton':
            res = optimize.newton(func, args=self.model.args, **kwargs)
        
        # multivariate methods are handled by optimize.root
        elif method in multivariate_methods:
            res = optimize.root(func, args=self.model.args, **kwargs)
        else:
            raise ValueError, 'Unrecognized method!'
    
        return res
        
    def set_eigen_star(self, eval_jacobian):
        """
        Computes the eigenvalues and left eigenvectors of the Jacobian matrix
        evaluated at a point in phase space. 
        
        Arguments:
            
            eval_jacobian: (array) Model jacobian evaluated at a point.
                        
        """
        # compute the eigenvalues and left eigenvectors 
        self.eigenvalues, self.eigenvectors = linalg.eig(eval_jacobian, 
                                                         left=True, 
                                                         right=False) 
        
    def check_local_stability(self, vec):
        """
        Checks whether local stability conditions are satisfied at a given 
        point. Stability conditions differ depending on whether or not the 
        timing of the model is continuous or discrete. In order for the 
        deterministic steady state to be locally stable we need to have:

            1. same number of stable eigenvalues as pre-determined 
               variables (i.e., state variables)
            2. same number of unstable eigenvalues as control (i.e., 
               jump variables).
        
        Definition of stable (unstable) eigenvalues differs slightly depending
        on whether the model is in continuous or discrete time. In continuous 
        time, a the real part of a stable (unstable) eigenvalue must have 
        negative (positive). In discrete time a stable (unstable) eigenvalue
        has modulus < (>) 1.
            
        Arguments:
            
            vec: (array) Point at which to conduct the local stability analysis.
          
        Returns: A list containing...

            eigenvalues:  The eigenvalues of the Jacobian matrix.
            eigenvectors: The eigenvectors of the Jacobian matrix.  

        """
        # evaluate the Jacobian
        eval_jacobian = self.model.evaluate_jacobian(vec)
        
        # compute eigenvalues and eigenvectors
        self.set_eigen_star(eval_jacobian)
