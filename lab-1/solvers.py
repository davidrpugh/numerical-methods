#import copy
import numpy as np
from numpy import polynomial 
from scipy import integrate, interpolate, optimize

class IVP(object):
    """Base class for solving initial value problems (IVPs)."""
    
    def __init__(self, f, jac, args=None):
        """
        Initializes an IVP object with the following attributes:
                    
            f:    (callable) Function returning the RHS of a system of ODEs. 
                                
            jac:  (callable) Function returning the model's jacobian matrix of 
                  partial derivatives. Must take the same arguments as f.
                         
            args: (tuple) Tuple of extra arguments which should be passed to the
                  functions f and jac. Default is None.
            
        """
        # initialize model attributes
        self.f    = f
        self.jac  = jac
        self.args = args

        # create and instance of the scipy.integrate.ode class      
        self.ode  = integrate.ode(f, jac)
            
    def integrate(self, t0, y0, h=1.0, T=None, g=None, tol=None, 
                  integrator='dopri5', step=False, relax=False, **kwargs):
        """
        Generates solution trajectories of the model given some initial 
        conditions.
        
        Arguments:
                
            t0:         (float) Initial condition for the independent variable.
            
            y0:         (float) Initial condition for the dependent variable. 
                                            
            h:          (float) Step-size for computing the solution.
            
            T:          (int) Length of desired trajectory.
            
            g:          (callable) Function of the form g(t, vec, f_args) that
                        provides stopping conditions for the integration. 
                        If specified, user must also specify a stopping 
                        tolerance, tol. Default is None. 
            
            tol:        (float) Stopping tolerance. On required if g is given.
                        Default is None. 
                        
            integrator: (str) Must be one of:
                        
                        'forward_euler':    Basic implementation of Euler's 
                                            method with fixed step size. See
                                            Judd (1998), Chapter 10, pg 341 for
                                            more detail.
                        
                        'backward_euler':   Basic implementation of the 
                                            implicit Euler method with a
                                            fixed step size.  See Judd (1998),
                                            Chapter 10, pg. 343 for more detail.
                        
                        'trapezoidal_rule': Basic implementation of the 
                                            trapezoidal rule with a fixed step 
                                            size.  See Judd (1998), Chapter 10, 
                                            pg. 344 for more detail.
                        
                        'erk2':             Second-order explicit Runge-Kutta.
                        
                        'erk3':             Third-order explicit Runge-Kutta.
                        
                        'erk4':             Fourth-order explicit Runge-Kutta.
                        
                        'erk5':             Fifth-order explicit Runge-Kutta. 
                        
                        'vode':             Real-valued Variable-coefficient ODE
                                            equation solver, with fixed leading
                                            coefficient implementation. It 
                                            provides implicit Adams method (for 
                                            non-stiff problems) and a method 
                                            based on backward differentiation 
                                            formulas (BDF) (for stiff problems).
                        
                        'lsoda':            Real-valued Variable-coefficient ODE
                                            equation solver, with fixed leading
                                            coefficient implementation. It 
                                            provides automatic method switching 
                                            between implicit Adams method (for 
                                            non-stiff problems) and a method
                                            based on backward differentiation 
                                            formulas (BDF) (for stiff problems).
                        
                        'dopri5':           Embedded explicit Runge-Kutta method 
                                            with order 4(5). See Dormand and 
                                            Prince (1980) for details. 
                        'dop85':
                                 
                        See documentation for integrate.ode for more details and 
                        references for 'vode', 'lsoda', 'dopri5', and 'dop85', 
                        as well as the rest of the ODE solvers available via
                        ODEPACK.
            
            step:       (boolean) Allows access to internal steps for those 
                        solvers that use adaptive step size routines. Currently
                        only 'vode', 'zvode', and 'lsoda' allow support step. 
                        Default is False. 
                         
            relax:      (boolean) The following integrators support run_relax: 
                        'vode', 'zvode', 'lsoda'. Default is False. 
                     
            **kwargs:   (dict) Dictionary of integrator specific keyword args.
                
        Returns: 
                     
           solution: (array-like) Simulated solution trajectory.
                            
        """        
        # select the integrator
        self.ode.set_integrator(integrator, **kwargs)
        
        # pass the model parameters as additional args
        if self.args != None:
            self.ode.set_f_params(*self.args)
            self.ode.set_jac_params(*self.args)
        
        # set the initial condition
        self.ode.set_initial_value(y0, t0)
        
        # create a storage container for the trajectory
        solution = np.hstack((t0, y0)) 
        
        # generate a solution trajectory 
        while self.ode.successful(): 
            self.ode.integrate(self.ode.t + h, step, relax)
            current_step = np.hstack((self.ode.t, self.ode.y))
            solution = np.vstack((solution, current_step))
                 
            # check terminal conditions
            if g != None and g(self.ode.t, self.ode.y, *self.args) < tol:
                break
            
            elif T != None and h > 0 and self.ode.t >= T:
                break
                
            elif T != None and h < 0 and self.ode.t <= T:
                break
            
            else:
                pass     
                            
        return solution
    
    def interpolate(self, traj, ti, k=3, der=0, ext=0):
        """
        Parameteric B-spline interpolation in N-dimensions.
        
        Arguments:
            
            traj: (array-like) Solution trajectory providing the data points for
                  constructing the B-spline representation.
                  
            ti:   (array-like) Array of values for the independent variable at 
                  which to interpolate the value of the B-spline.
            
            k:    (int) Degree of the desired B-spline. Degree must satsify 
                  1 <= k <= 5. Default is k=3 for cubic B-spline interpolation.
                  
            der:  (int) The order of derivative of the spline to compute 
                  (must be less than or equal to k). Default is zero.
            
            ext: (int) Controls the value of returned elements for outside the
                 original knot sequence provided by traj. For extrapolation, set
                 ext=0; ext=1 returns zero; ext=2 raises a ValueError. Default 
                 is to perform extrapolation.
                 
        Returns:
            
            interp_traj: (array) The interpolated trajectory.   
        
        """        
        # array of parameter values
        u = traj[:,0]
        
        # build list of input arrays
        n = traj.shape[1]
        x = [traj[:,i] for i in range(1, n)]
        
        # construct the B-spline representation (s=0 forces interpolation!)
        tck, t = interpolate.splprep(x, u=u, k=k, s=0)
        
        # evaluate the B-spline (returns a list)
        out = interpolate.splev(ti, tck, der, ext)
        
        # convert to a 2D array
        interp_traj = np.hstack((ti[:,np.newaxis], np.array(out).T))
        
        return interp_traj
             
    def compare_trajectories(self, traj1, traj2):
        """
        Returns the absolute difference between two solution trajectories.
        
        Arguments:
            
            traj1: (array-like) (T,n+1) array containing a solution trajectory.
            traj2: (array-like) (T,n+1) array containing a solution trajectory.
            
        Returns:
            
            abs_diff: (array-like) (T,n) array of the element-wise absolute 
                      difference between traj1 and traj2.
        """
        # element-wise absolute difference
        abs_diff = np.abs(traj1[:,1:] - traj2[:,1:])
        
        return abs_diff
        
    def get_L2_errors(self, traj1, traj2):
        """
        Computes a measure of the average difference between two trajectories
        using the L^2 norm.
        
        Arguments:
            
            traj1: (array-like) (T,n+1) array containing a solution trajectory.
            traj2: (array-like) (T,n+1) array containing a solution trajectory.
            
        Returns:
            
            L2_error: (float) Average difference between two trajectories.
                      
        """
        # L2 norm is the Euclidean distance
        L2_error = np.sum(self.compare_trajectories(traj1, traj2)**2)**0.5
        
        return L2_error
        
    def get_maximal_errors(self, traj1, traj2):
        """
        Computes a measure of the average difference between two trajectories
        using the L^oo norm.
        
        Arguments:
            
            traj1: (array-like) (T,n+1) array containing a solution trajectory.
            traj2: (array-like) (T,n+1) array containing a solution trajectory.
            
        Returns:
            
            maximal_error: (float) Maximal difference between two trajectories.
                      
        """
        # L^oo norm is the maximal distance
        maximal_error = np.max(self.compare_trajectories(traj1, traj2))
        
        return maximal_error   
                                    
class BVP(IVP):
    """Base class for solving boundary value problems (BVPs)."""

    def solve_generic_shooting(self, t0, y0, g0, R, R_args, method='hybr', 
                               options=None, h=1.0, T=None, tol=None, 
                               integrator='dopri5', **kwargs):
        """
        Solves a two-point BVP using a generic 'forward shooting' algorithm (see
        Judd (1998) p. 351 for details).

        Arguments:
            
            t0:          (float) Intitial condition for the independent 
                         variable.
                         
            y0:          (array-like) Given initial conditions for n=1,...n'
                         components of the solution, y(t).
            
            g0:          (array-like) Guess at the correct value for the 
                         remaining n - n' initial conditions for which only
                         terminal conditions are supplied.
                         
            R:           (callable) Function defining the boundary conditions.
            
            R_args:      (tuple) Additional arguments passed to the function R.
                         
            method:      (str) Valid method for solving non-linear equations.
            
            options:     (dict) Dictionary of method specific keyword arguments
                         passed to the nonlinear equation solver. 
            
            h:           (float) Step-size for computing the solution trajectory.
            
            tol:         (float) Convergence tolerance for solution trajectory. 
                         Due to accumulation of numerical error in solving the 
                         IVP for some k0 and c0, it may be necessary to choose a
                         relatively loose stopping criterion.
                         
            integrator:  (str) Must be a valid integrator class. See docstring 
                         of the integrate method for a complete listing of valid
                         integrators. 
                     
            **kwargs:    (dict) Dictionary of method specific keyword arguments
                         passed to the nonlinear equation solver. 
                        
        Returns: A list containing...
                     
            res       (object) Result object returned by scipy.optimize.root.
            
            solution: (array-like or None) Simulated solution trajectory 
                      approximating the solution to the BVP. If no suitable 
                      initial condition (i.e., one consistent with the boundary 
                      conditions) was found then this element will be None. 
               
        """      
        # first need to solve the non-linear equation for the optimal g0
        res = optimize.root(R, g0, args=R_args, method=method, options=options)

        if res.success == True:
            init = np.hstack((y0, res.x))
            solution = self.integrate(t0, init, h, T, integrator=integrator, 
                                      **kwargs)
            return [res, solution]
        else:
            return [res, None]
        
class bvpcol(IVP):
    """
    Base class for solving two-point boundary value problems using orthogonal 
    collocation methods.
    
    """
    
    def __init__(self, f, jac, a, b, bc, args=None):
        """
        Initializes a bvpcol object with the following attributes:
            
           f:    (callable) Function returning the RHS of a system of ODEs. 
                                
           jac:  (callable) Function returning the model's jacobian matrix of 
                  partial derivatives. Must take the same arguments as f.
                       
           a:    (float)
               
           b:    (float)
           
           bc:   (callable) Function of the form g(ya, yb, args) defining the 
                 boundary conditions.
           
           args: (tuple) Tuple of extra arguments which should be passed to the
                 functions f, jac, and bc. Default is None.
                 
        """
        super(bvpcol, self).__init__(f, jac, args)  
        self.a          = a
        self.b          = b
        self.bc         = bc
        
        self.yhat       = None
        self.yhat_prime = None
                      
    def set_yhat(self, coefs, kind='chebyshev'):
        """Modifies the yhat attribute."""
        if kind == 'polynomial':
            self.yhat, self.yhat_prime = self.__get_poly_yhat(coefs)
            
        elif kind == 'chebyshev':
            self.yhat, self.yhat_prime = self.__get_cheb_yhat(coefs)
            
        elif kind == 'legendre':
            self.yhat, self.yhat_prime = self.__get_lege_yhat(coefs)
            
        elif kind == 'laguerre':
            self.yhat, self.yhat_prime = self.__get_lagu_yhat(coefs)
            
        elif kind == 'hermite':
            self.yhat, self.yhat_prime = self.__get_herm_yhat(coefs)
            
        else:
            raise ValueError
   
    def __get_poly_yhat(self, coefs):
        """Constructs a standard polynomial approximation to the solution."""
        # the solution y is in R^n
        n = coefs.shape[0]
        
        # create the approximation yhat 
        domain = [self.a, self.b]
        yhat = [polynomial.Polynomial(coefs[i], domain) for i in range(n)]
        
        # create the approximation of yhat_prime
        yhat_prime = [yhat[i].deriv() for i in range(n)]
        
        return [yhat, yhat_prime]
         
    def __get_cheb_yhat(self, coefs):
        """Constructs a Chebyshev polynomial approximation to the solution."""
        # the solution y is in R^n
        n = coefs.shape[0]
        
        # create the approximation yhat 
        domain = [self.a, self.b]
        yhat   = [polynomial.Chebyshev(coefs[i], domain) for i in range(n)]
        
        # create the approximation of yhat_prime
        yhat_prime = [yhat[i].deriv() for i in range(n)]
        
        return [yhat, yhat_prime]
    
    def __get_lege_yhat(self, coefs):
        """Constructs a Legendre polynomial approximation to the solution."""
        # the solution y is in R^n
        n = coefs.shape[0]
        
        # create the approximation yhat 
        domain = [self.a, self.b]
        yhat = [polynomial.Legendre(coefs[i], domain) for i in range(n)]
        
        # create the approximation of yhat_prime
        yhat_prime = [yhat[i].deriv() for i in range(n)]
        
        return [yhat, yhat_prime]
    
    def __get_lagu_yhat(self, coefs):
        """Constructs a Laguerre polynomial approximation to the solution."""
        # the solution y is in R^n
        n = coefs.shape[0]
        
        # create the approximation yhat 
        domain = [self.a, self.b]
        yhat = [polynomial.Laguerre(coefs[i], domain) for i in range(n)]
        
        # create the approximation of yhat_prime
        yhat_prime = [yhat[i].deriv() for i in range(n)]
        
        return [yhat, yhat_prime]
          
    def __get_herm_yhat(self, coefs):
        """Constructs a Hermite polynomial approximation to the solution."""
        # the solution y is in R^n
        n = coefs.shape[0]
        
        # create the approximation yhat 
        domain = [self.a, self.b]
        yhat = [polynomial.Hermite(coefs[i], domain) for i in range(n)]
        
        # create the approximation of yhat_prime
        yhat_prime = [yhat[i].deriv() for i in range(n)]
        
        return [yhat, yhat_prime]
        
    def residual_function(self, t):
        """The residual function. Problem is here!"""
        # the solution is in R^n
        n = len(self.yhat)
         
        # evaluate the approximating functions
        y       = np.array([self.yhat[i](t) for i in range(n)])
        y_prime = np.array([self.yhat_prime[i](t) for i in range(n)])
        
        # compute the residual as (n, m) array
        residual = y_prime - self.f(t, y) # need to deal with case of args!
        
        return residual.flatten()
        
    def collocation_system(self, coefs, grid):
        """System of non-linear equations for collocation method."""
        # update the approximating functions coefs
        n = len(self.yhat)
        m = grid.size
        coefs = coefs.reshape((n, m + 1))
        
        for i in range(n):
            self.yhat[i].coefs       = coefs[i] 
            self.yhat_prime[i].coefs = coefs[i]
        
        # evaluate the residual functions
        ya = np.array([self.yhat[i](self.a) for i in range(n)])
        yb = np.array([self.yhat[i](self.b) for i in range(n)])
        
        out = np.hstack((self.bc(ya, yb), 
                         self.residual_function(grid)))
        return out
    
    def __get_collocation_nodes(self, m, kind):
        """Computes the collocation nodes."""        
        # basis coefs
        coefs = np.zeros(m + 1)
        coefs[-1] = 1
        
        if kind == 'polynomial':
            nodes  = polynomial.Polynomial(coefs).roots()
            
        elif kind == 'chebyshev':
            domain = [self.a, self.b]
            nodes  = polynomial.Chebyshev(coefs, domain).roots()
            
        elif kind == 'legendre':
            domain = [self.a, self.b]
            nodes  = polynomial.Legendre(coefs, domain).roots()
            
        elif kind == 'laguerre':
            domain = [self.a, self.b]
            nodes  = polynomial.Laguerre(coefs, domain).roots()
            
        elif kind == 'hermite':
            domain = [self.a, self.b]
            nodes  = polynomial.Hermite(coefs, domain).roots()
        
        else:
            raise ValueError

        return nodes
        
    def solve(self, coefs, kind='chebyshev', method='hybr', solver_opts=None):
        """
        Collocation method for solving 2PBVP.
        
        Arguments:
            
            coefs:  (array) An (n, m + 1) array of coefficients where n is the
                    dimension of the solution function and m is the desired 
                    degree of the approximating polynomial.
                   
            kind:   (str) Family of orthogonal polynomials to use as basis 
                    functions. Must be one of 'polynomial', 'chebyshev',
                    'legendre', 'laguerre', or 'hermit'. Default is 'chebyshev'.
                   
            method: (str) Method for solving the system of algebraic equations.
                    Must be one of 'hybr', 'lm', 'broyden1', 'broyden2', 
                    'anderson', 'linearmixing', 'diagbroyden', 'excitingmixing',
                    'krylov'. See the scipy.optimize.root documentation for
                    solver specific details.
            
        """
        # create the approximating functions
        self.set_yhat(coefs, kind)
        
        # compute the collocation nodes
        n = coefs.shape[0]
        m = coefs.shape[1] - 1
        nodes = self.__get_collocation_nodes(m, kind)
        
        res = optimize.root(self.collocation_system, coefs.flatten(), 
                            args=(nodes,), method=method, options=solver_opts) 
        
        if res.success == False:
            print 'Booo!!!'
            return res
            
        else:
            print 'Yeah!!!'
            solution = res.x.reshape((n, m + 1))
            self.set_yhat(solution, kind)        
