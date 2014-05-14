from __future__ import division

import pandas as pd
import numpy as np
from scipy import integrate, interpolate, optimize

import matplotlib as mpl
import matplotlib.pyplot as plt

import solvers, integrators
import steady_states
import pwt
          
class SolowModel(solvers.IVP):
    """Base class for the Solow (1956) model of growth."""
    pwt_data = pwt.load_pwt_data()
        
    def __init__(self, output, mpk, k_dot, jacobian, params=None):
        """
        Initializes a SolowModel object with the following attributes:
            
            output:   (callable) Output (per person/effective person). Should be
                      of the form
                      
                          f(t, k, params)
                          
                      where the independent variable, t, is time; k, is capital
                      per effective person; and params is a dictionary of model
                      parameters.       
            
            mpk:      (callable) Marginal product of capital (per person/
                      effective person). Should be of the form
                      
                          f'(t, k, params)
                          
                      where the independent variable, t, is time; k, is capital
                      per effective person; and params is a dictionary of model
                      parameters.
                         
            k_dot:    (callable) Equation of motion for capital (per person/
                      effective person). Should be of the form 
                     
                         k_dot(t, k, params)
                        
                      where the independent variable, t, is time; k, is capital
                      per effective person; and params is a dictionary of model
                      parameters.
            
            jacobian: (callable) Returns the Jacobian matrix of partial 
                      derivatives for F. Should be of the form
                      
                          jacobian(t, vec, params)
                       
                      where the independent variable, t, is time; k, is capital
                      per effective person; and params is a dictionary of model
                      parameters.
                                   
            params:  (dict) Dictionary of model parameters. Default is None.
            
        """
        # initialize model attributes
        self.output                   = output
        self.mpk                      = mpk
        self.k_dot                    = k_dot
        self.jacobian                 = jacobian
        self.iso3_code                = None
        self.data                     = None
        self.dep_rates                = None
        
        # initialize an empty SteadyState object
        self.steady_state = steady_states.SteadyState(self) 
        
        # initialize the model as an IVP
        super(SolowModel, self).__init__(self.k_dot, self.jacobian, (params,))            
        self.params = self.args[0]
                
    def update_model_parameters(self, new_params):
        """Updates the model's parameter dictionary."""
        self.args = (new_params.copy(),) + self.args[1:]
        self.params = self.args[0]
       
    def linearized_k_dot(self, k):
        """
        Linearized version of the equation of motion for capital (per person/
        effective person).
        
        Arguments:
            
            k: (array-like) Capital (per person/effective person).
            
        Returns:
            
            linear_k_dot: (array) Linearized equation of motion for capital
                          (per person/effective person).
        
        """
        # evaluate the jacobian at steady state 
        k_star = self.steady_state.values['k_star']    
        lmbda = self.jacobian(0, k_star, self.params)
    
        # linearize!
        linear_k_dot = lmbda * (k - k_star)
    
        return linear_k_dot  
                 
    def linearized_solution(self, k0, t):
        """
        Linear approximation of the time-path of capital per (person/
        effective person).
        
        Arguments:
            
            k0: (float) Initial condition for capital (per person/effective 
                person).
                
            t:  (array) (T,) array of time points at which to return a solution.
        
        Returns:
            
            linear_traj: (array) (T,2) array representing a linearized time-path
                         for capital (per person/effective person).
                
        """
        # get steady state value of capital
        k_star = self.steady_state.values['k_star']
        
        # speed of convergence
        lmbda = self.jacobian(0, k_star, self.params)
        k_linear = k_star + np.exp(lmbda * t) * (k0 - k_star)
        
        # construct the linearized trajectory
        linear_traj = np.hstack((t[:, np.newaxis], k_linear[:, np.newaxis]))
        
        return linear_traj
    
    def get_capitals_share(self, k, params):
        """
        Computes capital's share of income alpha_k(k) for the Solow model. 
        
        """
        alpha_k = (k * self.mpk(0, k, params)) / self.output(0, k, params)
    
        return alpha_k
                
    def get_impulse_response(self, param, shock, T=100, year=2013, 
                             kind='efficiency_units', reset=True):
        """
        Generates an impulse response function for k(t) following a shock to one
        of the model parameters.
        
        Arguments:
            
            param: (string) Model parameter
            
            shock: (float) Multiplicative shock to the parameter. Values < 1 
                   correspond to a reduction in the current value of the 
                   parameter; values > 1 correspond to increasing the current 
                   value of the parameter.
                   
            T:     (float) Length of the impulse response. Default is 100.
            
            kind:  (str) Whether you want impulse response functions in 'levels',
                   'per_capita', or 'efficiency_units'. Default is for irfs to
                   be in 'efficiency_units'.
                   
            reset: (boolean) Whether or not to reset the original parameters to
                   their pre-shock values. Default is True.
            
        Returns:
            
            irf: (array-like) Impulse response function.
            
        """
        # copy the original params
        orig_params = self.params.copy()
        
        # economy is initial in steady state
        k0 = self.steady_state.values['k_star']
        y0 = self.output(0, k0, orig_params)
        c0 = (1 - orig_params['s']) * y0
        
        # initial padding should be such that shock occurs in year
        timedelta = pd.datetime(year, 1, 1) - self.data.index[0]
        N = timedelta.days // 365
        time_padding = np.arange(0, N, 1.0)
        
        # transform irfs into per capita or levels, depending
        if kind == 'per_capita':
            A0 = self.params['A0']
            g  = self.params['g']
            factor = A0 * np.exp(g * time_padding)
            
        elif kind == 'levels':
            A0 = self.params['A0']
            g  = self.params['g']
            L0 = self.params['L0']
            n  = self.params['n']
            factor = A0 * L0 * np.exp((n + g) * time_padding)
            
        elif kind == 'efficiency_units':
            factor = np.ones(N)
        
        else:
            raise ValueError
            
        # start with N periods of steady state values
        padding_k = np.repeat(k0, N)
        padding   = np.hstack((time_padding[:,np.newaxis], 
                              (factor * padding_k)[:,np.newaxis]))
        # padding for y
        padding_y = np.repeat(y0, N)
        padding   = np.hstack((padding, (factor * padding_y)[:,np.newaxis]))
        
        # padding for c
        padding_c = np.repeat(c0, N)
        padding   = np.hstack((padding, (factor * padding_c)[:,np.newaxis]))
            
        # shock the parameter
        self.params[param] = shock * self.params[param]
        self.update_model_parameters(self.params)
        
        # compute the new steady state values
        self.steady_state.set_values()
        
        # generate post-shock trajectory
        irf = self.integrate(1.0, k0, 1.0, T, integrator='dopri')     
        
        # transform irfs into per capita or levels, depending
        if kind == 'per_capita':
            g  = self.params['g']
            factor = factor[-1] * np.exp(g * irf[:,0])
            
        elif kind == 'levels':
            g  = self.params['g']
            n  = self.params['n']
            factor = factor[-1] * np.exp((n + g) * irf[:,0])
            
        elif kind == 'efficiency_units':
            factor = np.ones(T)
        
        else:
            raise ValueError
        
        # compute the irf for y
        irf_y = self.output(irf[:,0], irf[:,1], self.params)
        irf = np.hstack((irf, (factor * irf_y)[:,np.newaxis]))
        
        # compute the irf for c
        irf_c = (1 - self.params['s']) * irf_y
        irf = np.hstack((irf, (factor * irf_c)[:,np.newaxis]))
        
        # compute the irf for k (after computing irfs for y and c!)
        irf[:,1] = (factor * irf[:,1])
        
        # add the padding
        irf = np.vstack((padding, irf))
        
        # shift the time index forward
        irf[:N, 0] += self.data.index[0].year 
        irf[N:, 0] += year - 1
        
        # reset the original params and recompute steady state?
        if reset == True:
            self.update_model_parameters(orig_params)
            self.steady_state.set_values()
        
        return irf      
                                 
    def plot_approximation_error(self, numeric_traj, analytic_traj, log=False):
        """
        Plots the numerical approximation error.
        
        Arguments:
            
            numeric_traj:  (array-like) (T,2) array containing a numeric 
                           trajectory.
            analytic_traj: (array-like) (T,2) array containing an analytic 
                           trajectory.
            log:           (boolean) Whether or not you wish to have a log-scale 
                           for the vertical axis. Default is False.
            
        Returns: A list containing...
        
            ax:   (object) Axes object representing the plot.
            line: (object) Line2D object representing the approximation error.
        
        """
        # create a new figure and subplot
        ax     = plt.subplot(111)
        
        # plot the approximation error
        approx_error = self.compare_trajectories(numeric_traj, analytic_traj)
        line         = ax.plot(numeric_traj[1:,0], approx_error[1:])[0]

        # logarithmic scale for the y-axis?
        if log:
            ax.set_yscale('log')
            ax.set_ylabel(r'$|k_n - k(t_n)|$ (log scale)', family='serif', 
                          fontsize=15)
        else:
            ax.set_ylabel(r'$|k_n - k(t_n)|$', family='serif', fontsize=15)
        
        # set the x-axis label
        ax.set_xlabel('Time, $t$', family='serif', fontsize=15)
            
        # remove the right and top spines
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # hide the top and right ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # provide a title
        ax.set_title('Approximation error for $k(t)$', family='serif', 
                     fontsize=20)

        return [ax, line]
        
    def plot_phase_diagram(self, N=1000, linearize=False):
        """
        Generates a plot of the phase diagram for the Solow model.
        
        Arguments:
            
            N: (int, optional) Number of points to plot. Default is 1000. 
            
        Returns: A list containing...
        
            ax:   (object) Axes object representing the plot.
            line: (object) Line2D object representing the k-dot locus.
        
        """
        # create a new figure and subplot
        ax     = plt.subplot(111)

        # use value of k_star to anchor the plot
        k_star = self.steady_state.values['k_star']
        grid   = np.linspace(0, 2 * k_star, N)

        # plot the evolution of capital
        data   = self.f(0, grid, *self.args)
        line   = ax.plot(grid, data)

        # add linearize equation of motion
        if linearize == True:
            ax.plot(grid, self.linearized_k_dot(grid))
            
        # adjust the ylims
        ax.set_ylim(-np.max(np.abs(data)), np.max(np.abs(data)))

        # remove the right and top spines
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # hide the top and right ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # center the bottom spine on the data and demarcate the steady state
        ax.spines['bottom'].set_position(('data', 0))
        ax.set_xticks([k_star])
        ax.set_xticklabels([r'$k^*$'])

        # label the vertical axis
        ax.set_ylabel(r'$\dot{k}$', rotation='horizontal', fontsize=15, 
                      family='serif')

        # provide a title
        ax.set_title('Phase diagram for $k$ in the Solow model', fontsize=20, 
                     family='serif')

        return [ax, line]
        
    def plot_solow_diagram(self, gridmax, N=1000, param=None, shock=None, 
                           reset=True):
        """
        Generates the classic Solow diagram.
        
        Arguments:
            
            gridmax: (float) Maximum value for capital per effective person.
            
            N:       (int, optional) Number of points to plot. Default is 1000.
             
            param:   (string) Model parameter.
            
            shock:   (float) Multiplicative shock to the parameter. Values < 1 
                     correspond to a reducing the current value of a parameter; 
                     values > 1 correspond to increasing the current value of 
                     the parameter.
                     
            reset:   (boolean) Whether or not to reset the original parameters
                     to their pre-shock values. Default is True.
        
        Returns: A list containing...
        
            ax:          (object) Axes object representing the plot.
            output:      (object) Line2D object representing output.
            act_inv:     (object) Line2D object representing actual investment.
            br_even_inv: (object) Line2D object representing break-even
                         investment.
        
        """
        # extract params
        n     = self.params['n']
        g     = self.params['g']
        s     = self.params['s']
        delta = self.params['delta']
        
        # create a new figure and subplot
        ax    = plt.subplot(111)

        # grid of values for capital per effective person
        grid  = np.linspace(0, gridmax, N)
          
        # plot output, actual and break even investment             
        output = ax.plot(grid, self.output(0, grid, self.params), 'r', 
                         label='$y$')[0]
        act_inv = ax.plot(grid, s * self.output(0, grid, self.params), 'g', 
                          label='$i_{act}$')[0]
        br_even_inv = ax.plot(grid, (n + g + delta) * grid, 'b', 
                              label='$i_{br}$')[0]
           
        # remove the right and top spines
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # hide the top and right ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # axes, labels, title, legend, etc
        ax.set_xlabel('Capital per effective person, $k$', fontsize=15, 
                      family='serif')
        ax.set_xlim(0, gridmax)
       
        # handles parameter shocks   
        if param != None and shock !=None:
            
            # copy the original params
            orig_params = self.params.copy()
            
            # shock the parameter
            self.params[param] = shock * self.params[param]
            self.update_model_parameters(self.params)
            
            # compute the new steady state values
            self.steady_state.set_values()
            
            # extract possibly new params
            n     = self.params['n']
            g     = self.params['g']
            s     = self.params['s']
            delta = self.params['delta']
            
            # plot new output, actual and break even investment             
            new_output = ax.plot(grid, self.output(0, grid, self.params),'r')[0]
            new_act_inv = ax.plot(grid, s * self.output(0, grid, self.params),
                                  'g')[0]
            new_br_even_inv = ax.plot(grid, (n + g + delta) * grid, 'b')[0]
            
            # plot formatting depends on parameter being shocked
            if param in ['alpha', 'delta', 'sigma']:
                param = '\\' + param  # necessary for pretty latex printing
                
            if param in ['n', 'g', '\\delta']:
                br_even_inv.set_alpha(0.5)
                br_even_inv.set_label(r'$i_{br, old}$')
                new_br_even_inv.set_label(r'$i_{br, new}$')
                ax.set_title('Changing $%s$ shifts break-even investment!' 
                             %param, fontsize=20, family='serif')
                ax.legend(loc='best', frameon=False)
                
            elif param == 's':
                act_inv.set_alpha(0.5)
                act_inv.set_label(r'$i_{act, old}$')
                new_act_inv.set_label(r'$i_{act, new}$')
                ax.set_title('Changing $%s$ shifts actual investment!' %param, 
                             fontsize=20, family='serif')
                ax.legend(loc='best', frameon=False) 
            
            elif param in ['\\alpha', '\\sigma']:
                output.set_alpha(0.5)
                output.set_label(r'$y_{old}$')
                new_output.set_label(r'$y_{new}$')
                act_inv.set_alpha(0.5)
                act_inv.set_label(r'$i_{act, old}$')
                new_act_inv.set_label(r'$i_{act, new}$')
                ax.set_title('Changing $%s$ shifts output and actual investment!' 
                             %param, fontsize=20, family='serif')
                ax.legend(loc='best', frameon=False)
                
            else:
                raise ValueError
            
            # reset the original params and recompute steady state?
            if reset == True:
                self.update_model_parameters(orig_params)
                self.steady_state.set_values()
            
            out = [ax, output, act_inv, br_even_inv, new_output, new_act_inv, 
                   new_br_even_inv]
                   
        else:
            ax.set_title('Classic Solow Diagram\n' + 
                         '$s=%g, n=%g, g=%g, \delta=%g$' %(s,n,g,delta), 
                         fontsize=20, family='serif')
            ax.legend(loc='best', frameon=False) 
            
            out = [ax, output, act_inv, br_even_inv]
               
        return out
        
    def plot_trajectory(self, traj, color='b', ax=None, **kwargs):
        """
        Plots the time path of the economy.
        
        Arguments:
            
            traj:     (array) Array representing a time path of the economy. 
                      If `kind` is set to `time_series', then the first column 
                      of traj should be the time index against which any 
                      additional columns should be plotted. 
            
            color:    (varies) Valid matplotlib color specification.
            
            ax:       (object) AxesSubplot object on which the trajectory 
                      should be plotted.
                               
            **kwargs: (dict) Dictionary of optional keyword arguments to pass
                      to ax.plot method.
                   
        """ 
        if ax == None:
            ax = plt.subplot(111)
                
        ax.plot(traj[:,0], traj[:,1], color=color, **kwargs)
        
    def plot_impulse_response(self, variables, param, shock, T, year=2013,
                              color='b', kind='efficiency_units', log=False, 
                              reset=True, **fig_kw):
        """
        Plots an impulse response function.
        
        Arguments:
            
            variables: (list) List of variables whose impulse response functions
                       you wish to plot. Alternatively, you can plot irfs for 
                       all variables by setting variables = 'all'.
                    
            param:     (string) Model parameter you wish to shock.
            
            shock:     (float) Multiplicative shock to the parameter. Values < 1 
                       correspond to a reduction in the current value of the 
                       parameter; values > 1 correspond to an increase in the 
                       current value of the parameter.
                   
            T:         (float) Length of the impulse response. Default is 40.
            
            year:      (int) Year in which you want the shock to take place.
                       Default is 2013.
                       
            kind:      (str) Whether you want impulse response functions in 
                       'levels', 'per_capita', or 'efficiency_units'. Default is
                       for irfs to be in 'efficiency_units'.
                   
            log:       (boolean) Whether or not to have logarithmic scales on
                       the vertical axes. Default is False.
                   
            reset:     (boolean) Whether or not to reset the original parameters
                       to their pre-shock values. Default is True.
            
        Returns: A list containing...
        
            fig:  (object) Instance of the matplotlib Figure class
            axes: (list) List of matplotlib AxesSubplot instances.
        
        """
        # first need to generate and irf
        irf = self.get_impulse_response(param, shock, T, year, kind, reset)
        
        # create mapping from variables to column indices
        irf_dict = {'k':irf[:,[0,1]], 'y':irf[:,[0,2]], 'c':irf[:,[0,3]]}
        
        # necessary for pretty latex printing
        if param in ['alpha', 'delta', 'sigma', 'theta', 'rho']:
            param = '\\' + param
            
        # title depends on whether shock was positive or negative
        if shock > 1.0:
            tit = 'Impulse response following + shock to $%s$' % param
        elif shock < 1.0:
            tit = 'Impulse response following - shock to $%s$' % param
        else:
            tit = 'Impulse response following NO shock to $%s$' % param

        if variables == 'all':
            variables = irf_dict.keys()
            
        fig, axes = plt.subplots(len(variables), 1, squeeze=False, **fig_kw)
          
        for i, var in enumerate(variables): 
                
            # extract the time series
            traj = irf_dict[var]
                
            # plot the irf
            self.plot_trajectory(traj, color, axes[i,0])
                
            # adjust axis limits
            axes[i,0].set_ylim(0.95 * traj[:,1].min(), 1.05 * traj[:,1].max())
            axes[i,0].set_xlim(year - 10, year + T)
                
            # y axis labels depend on kind of irfs
            if kind == 'per_capita':
                  
                ti = traj[:,0] - self.data.index[0].year
                gr = self.params['g']
                    
                axes[i,0].plot(traj[:,0], traj[0,1] * np.exp(gr * ti), 'k--')
                axes[i,0].set_ylabel(r'$\frac{%s}{L}(t)$' %var.upper(), 
                                     rotation='horizontal', fontsize=15, 
                                     family='serif')
                                               
            elif kind == 'levels':
                ti = traj[:,0] - self.data.index[0].year
                gr = self.params['n'] + self.params['g']
                    
                axes[i,0].plot(traj[:,0], traj[0,1] * np.exp(gr * ti), 'k--')
                axes[i,0].set_ylabel('$%s(t)$' %var.upper(), 
                                     rotation='horizontal', fontsize=15, 
                                     family='serif')
                                               
            else:
                axes[i,0].set_ylabel('$%s(t)$' %var, rotation='horizontal', 
                                     fontsize=15, family='serif')
                                       
            # adjust location of y-axis label
            axes[i,0].yaxis.set_label_coords(-0.1, 0.45)
    
            # log the y-scale for the plots
            if log == True:
                axes[i,0].set_yscale('log')
                    
        axes[0,0].set_title(tit, family='serif', fontsize=20)
        axes[-1,0].set_xlabel('Year, $t$,', fontsize=15, family='serif')
        
        return [fig, axes]
    
def calibrate_cobb_douglas(model, iso3_code, bounds=None):
    """
    Calibrates a Solow model with Cobb-Douglas production using data from the 
    Penn World Tables (PWT).

    Arguments:
        
        model:     (object) An instance of the SolowModel class.
            
        iso3_code: (str) A valid ISO3 country code.
        
        bounds:    (tuple) Start and end years for the subset of the PWT data
                   to use for calibration.
                    
    """
    # modify the country attribute
    model.iso3_code = iso3_code
        
    # get the PWT data for the iso_code
    model.data      = model.pwt_data.major_xs(iso3_code)
    
    # set bounds
    if bounds == None:
        start = model.data.index[0]
        end   = model.data.index[-1]
    else:
        start = bounds[0]
        end   = bounds[1]
        
    ##### estimate capital's share of income/output ####
    alpha = (1 - model.data.labsh.loc[start:end]).mean()
            
    ##### estimate the fraction of output saved #####
    s = model.data.csh_i.loc[start:end].mean()
        
    ##### estimate the labor force growth rate #####
        
    # regress log employed persons on linear time trend
    N     = model.data.index.size
    trend = pd.Series(np.linspace(0, N - 1, N), index=model.data.index)
    res   = pd.ols(y=np.log(model.data.emp.loc[start:end]), 
                   x=trend.loc[start:end])
    n     = res.beta[0]
    L0    = np.exp(res.beta[1])
    
    ##### estimate the technology growth rate #####
        
    # adjust measure of TFP 
    model.data['atfpna'] = (model.data.rtfpna**(1 / model.data.labsh) * 
                            model.data.hc)
        
    # regress log TFP on linear time trend
    res   = pd.ols(y=np.log(model.data.atfpna.loc[start:end]), 
                   x=trend.loc[start:end])
    g     = res.beta[0]
    A0    = np.exp(res.beta[1])
                           
    ##### estimate the depreciation rate for total capital #####
        
    # use average depreciation rate for total capital
    delta = model.data.delta_k.loc[start:end].mean()   
         
    # create a dictionary of model parameters
    params = {'s':s, 'alpha':alpha, 'delta':delta, 'n':n, 'L0':L0, 'g':g,
              'A0':A0, 'sigma':1.0}
        
    # update the model's parameters
    model.update_model_parameters(params)
                    
    # compute new steady state values
    model.steady_state.set_values()

def calibrate_ces(model, iso3_code, x0, base=2005, bounds=None, method='COBYLA', 
                  tol=1e-9, **kwargs):
    """
    Calibrates a Solow model with constant elasticity of substition (CES)
    production using data from the Penn World Tables (PWT).

    Arguments:
        
        model:     (object) An instance of the SolowModel class.
            
        iso3_code: (str) A valid ISO3 country code.        
                   
        x0:        (array-like) Initial guess of optimal alpha and sigma.
        
        method:    (str) Method used to estimate alpha and sigma via non-linear 
                   least-squares. Default is 'COBYLA'.
                   
        tol:       (float) Convergence tolerance for non-linear least squares 
                   problem. Default is 1e-9.
                   
        kwargs:    (dict) Dictionary of additional keyword options to pass to the 
                   non-linear optimization routine.
    
    Returns: A list containing...
    
        res1: (object) OLS regression result used in estimating population growth rate.
        res2: (object) OLS regression result used in estimating the growth rate of
                       productivity.
        res3: (object) NNLS estimation result used to estimate alpha and sigma.
                
    """
    # modify the country attribute
    model.iso3_code = iso3_code
        
    # get the PWT data for the iso_code
    model.data      = model.pwt_data.major_xs(iso3_code)
    
    # set bounds for calibration data
    if bounds == None:
        start = model.data.index[0]
        end   = model.data.index[-1]
    else:
        start = bounds[0]
        end   = bounds[1]
    
    # extract model data used in calibration
    K       = model.data.rkna.loc[start:end]
    L       = model.data.emp.loc[start:end]
    Y       = model.data.rgdpna.loc[start:end]
    I_share = model.data.csh_i.loc[start:end]
    K_dep   = model.data.delta_k.loc[start:end]
    
    # compute adjusted measure of TFP 
    model.data['atfpna'] = (model.data.rtfpna**(1 / model.data.labsh) * 
                            model.data.hc)
    
    ##### estimate the fraction of output saved #####
    
    # use average investment share of gdp
    s = I_share.mean()
        
    ##### estimate the labor force growth rate #####
        
    # regress log employed persons on linear time trend
    N     = model.data.index.size
    trend = pd.Series(np.linspace(0, N - 1, N), index=model.data.index)
    res1  = pd.ols(y=np.log(L), x=trend.loc[start:end])
    n     = res1.beta[0]
    L0    = np.exp(res1.beta[1])
    
    ##### estimate the technology growth rate #####
        
    # regress log TFP on linear time trend
    A    = model.data.atfpna.loc[start:end]
    res2 = pd.ols(y=np.log(A), x=trend.loc[start:end])
    g    = res2.beta[0]
    A0   = np.exp(res2.beta[1])
    
    ##### estimate the depreciation rate for total capital #####
        
    # use average depreciation rate for total capital
    delta = K_dep.mean()  
       
    ##### estimate alpha and sigma using NLLS #####
    
    def objective(ces_params):
        """
        Choose parameters of ces producting function, alpha and sigma, 
        to minimize this function.
        
        """
        # extract parameters
        alpha, sigma = ces_params
        rho = (sigma - 1) / sigma
        
        # define the objective function
        y = (Y / (A * L))
        k = (K / (A * L))
    
        if abs(rho) < 1e-6:
            residual = np.log(y) - alpha * np.log(k)
            tss = np.sum(residual**2)
            
        else:
            residual = np.log(y) - (1 / rho) * np.log(alpha * k**rho + (1 - alpha))
            tss = np.sum(residual**2)
            
        return tss
    
    def constraint1(ces_params):
        """Insures that steady-state capital per effective worker is finite."""
        # extract parameters
        alpha, sigma = ces_params
        rho = (sigma - 1) / sigma
        
        return min(((n + g + delta) / s)**rho, 1) - alpha 
        
    def constraint2(ces_params):
        """Non-negativiely constraint on sigma."""
        # extract parameters
        alpha, sigma = ces_params  
        return sigma 
    
    def constraint3(ces_params):
        """Non-negativiely constraint on alpha."""
        # extract parameters
        alpha, sigma = ces_params        
        return alpha
    
    # create a dictionary of constraints
    constraints = [{'type':'ineq', 'fun':constraint1},
                   {'type':'ineq', 'fun':constraint2},
                   {'type':'ineq', 'fun':constraint3}]
    
    # solve the non-linear least squares problem        
    res3 = optimize.minimize(objective, x0, method=method, tol=tol,
                             constraints=constraints, **kwargs)
    alpha, sigma = res3.x
        
    # create a dictionary of model parameters
    params = {'s':s, 'alpha':alpha, 'sigma':sigma, 'delta':delta, 'n':n, 
              'L0':L0, 'g':g, 'A0':A0}
              
    # update the model's parameters
    model.update_model_parameters(params)
                    
    # compute new steady state values
    model.steady_state.set_values()

    return [res1, res2, res3]
