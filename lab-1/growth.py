from __future__ import division

import pandas as pd
import numpy as np
from scipy import integrate, interpolate, optimize
#import scikits.bvp_solver as bvp_solver

import matplotlib as mpl
import matplotlib.pyplot as plt

import solvers, integrators
import steady_states
import pwt
          
class SolowModel(solvers.IVP):
    """Base class for the Solow (1956) model of growth."""
    pwt_data, pwt_dep_rates = pwt.load_pwt_data(deltas=True)
        
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
                
    def get_impulse_response(self, param, shock, T=40, kind='efficiency_units', 
                             reset=True):
        """
        Generates an impulse response function for k(t) following a shock to one
        of the model parameters.
        
        Arguments:
            
            param: (string) Model parameter
            
            shock: (float) Shock to the parameter. Values < 1 correspond to a
                   reducing the current value of the parameter; values > 1
                   correspond to increasing the current value of the parameter.
                   
            T:     (float) Length of the impulse response. Default is 40.
            
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
        
        # initial padding should be 2000 time 'periods'
        time_padding = np.arange(0, 2001, 1.0)
        
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
            factor = np.ones(2001)
        
        else:
            raise ValueError
            
        # start with 2001 periods of steady state values
        padding_k = np.repeat(k0, 2001)
        padding   = np.hstack((time_padding[:,np.newaxis], 
                              (factor * padding_k)[:,np.newaxis]))
        # padding for y
        padding_y = self.output(padding[:,0], padding_k, self.params)
        padding = np.hstack((padding, (factor * padding_y)[:,np.newaxis]))
        
        # padding for c
        padding_c = (1 - self.params['s']) * padding_y
        padding = np.hstack((padding, (factor * padding_c)[:,np.newaxis]))
            
        # shock the parameter
        self.params[param] = shock * self.params[param]
        self.update_model_parameters(self.params)
        
        # compute the new steady state values
        self.steady_state.set_values()
        
        # generate post-shock trajectory
        irf = self.integrate(1, k0, 1.0, T, integrator='dopri')     
        
        # transform irfs into per capita or levels, depending
        if kind == 'per_capita':
            #A0 = self.params['A0']
            g  = self.params['g']
            
            factor = factor[-1] * np.exp(g * irf[:,0])
            
        elif kind == 'levels':
            #A0 = self.params['A0']
            g  = self.params['g']
            #L0 = self.params['L0']
            n  = self.params['n']
            
            factor = factor[-1] * np.exp((n + g) * irf[:,0])
            
        elif kind == 'efficiency_units':
            factor = np.ones(T)
        
        else:
            raise ValueError
        
        # shift the time index forward 2000 periods 
        irf[:,0] += 2000
        
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
       
    def plot_impulse_response(self, param, shock, T, kind='efficiency_units', 
                              log=False, reset=True):
        """
        Plots an impulse response function.
        
        Arguments:
            
            param: (string) Model parameter.
            
            shock: (float) Shock to the parameter. Values < 1 correspond to a
                   reducing the current value of the parameter; values > 1
                   correspond to increasing the current value of the parameter.
                   
            T:     (float) Length of the impulse response. Default is 40.
            
            kind:  (str) Whether you want impulse response functions in 'levels',
                   'per_capita', or 'efficiency_units'. Default is for irfs to
                   be in 'efficiency_units'.
                   
            log:   (boolean) Whether or not to have logarithmic scales on the
                   vertical axes. Default is False.
                   
            reset: (boolean) Whether or not to reset the original parameters to
                   their pre-shock values. Default is True.
            
        Returns: A list containing...
            
            ax1:   (Axes) Axes object for the impulse response for capital.
            irf_k: (Line2D) Impulse response for capital.
            ax2:   (Axes) Axes object for the impulse response for output.
            irf_y: (Line2D) Impulse response for output.
            ax3:   (Axes) Axes object for the impulse response for consumption.
            irf_c: (Line2D) Impulse response for consumption.
            
        TODO: Better way to deal with logarithmic y-scale (i.e., Don't force a
        common scale!)
        
        """
        # first need to generate and irf
        irf = self.get_impulse_response(param, shock, T, kind, reset)
        
        # irf for capital
        ax1 = plt.subplot(311)
        irf_k = ax1.plot(irf[:,0], irf[:,1], 'g')[0]
 
        ax1.set_xlim(1950, 2000 + T)
        ax1.yaxis.set_label_coords(-0.055, 0.45)
        
        # necessary for pretty latex printing
        if param in ['alpha', 'delta', 'sigma']:
            param = '\\' + param
            
        # title depends on whether shock was positive or negative
        if shock > 1.0:
            tit = 'Impulse response following + shock to $%s$' % param
        elif shock < 1.0:
            tit = 'Impulse response following - shock to $%s$' % param
        else:
            tit = 'Impulse response following NO shock to $%s$' % param
            
        ax1.set_title(tit, family='serif', fontsize=20)

        # irf for output
        ax2 = plt.subplot(312)
        irf_y = ax2.plot(irf[:,0], irf[:,2], 'r')[0]

        ax2.set_xlim(1950, 2000 + T)
        ax2.yaxis.set_label_coords(-0.055, 0.45)

        # irf for consumption
        ax3 = plt.subplot(313)
        irf_c = ax3.plot(irf[:,0], irf[:,3], 'b')[0]

        ax3.set_xlim(1950, 2000 + T)
        ax3.set_xlabel('Time, $t$,', fontsize=15, family='serif')
        ax3.yaxis.set_label_coords(-0.055, 0.45)

        if kind == 'efficiency_units':
            ax1.set_ylabel('$k(t)$', rotation='horizontal', fontsize=15, 
                       family='serif')
            ax2.set_ylabel('$y(t)$', rotation='horizontal', fontsize=15, 
                           family='serif')
            ax3.set_ylabel('$c(t)$', rotation='horizontal', fontsize=15, 
                           family='serif')
                           
        elif kind == 'per_capita':
            ax1.set_ylim(1e0, 1e3)
            ax2.set_ylim(1e0, 1e3)
            ax3.set_ylim(1e0, 1e3)
            
            ax1.set_ylabel(r'$\frac{K}{L}(t)$', rotation='horizontal', 
                           fontsize=15, family='serif')
            ax2.set_ylabel(r'$\frac{Y}{L}(t)$', rotation='horizontal', 
                           fontsize=15, family='serif')
            ax3.set_ylabel(r'$\frac{C}{L}(t)$', rotation='horizontal', 
                           fontsize=15, family='serif')
                           
        elif kind == 'levels':
            ax1.set_ylim(1e0, 5e5)
            ax2.set_ylim(1e0, 5e5)
            ax3.set_ylim(1e0, 5e5)
            
            ax1.set_ylabel('$K(t)$', rotation='horizontal', fontsize=15, 
                       family='serif')
            ax2.set_ylabel('$Y(t)$', rotation='horizontal', fontsize=15, 
                           family='serif')
            ax3.set_ylabel('$C(t)$', rotation='horizontal', fontsize=15, 
                           family='serif')
                           
        # log the y-scale for the plots
        if log == True:
            ax1.set_yscale('log')
            ax2.set_yscale('log')
            ax3.set_yscale('log')
            
        return [ax1, irf_k, ax2, irf_y, ax3, irf_c]

def calibrate_cobb_douglas(model, iso3_code):
    """
    Calibrates a Solow model with Cobb-Douglas production using data from the 
    Penn World Tables (PWT).

    Arguments:
        
        model:     (object) An instance of the SolowModel class.
            
        iso3_code: (str) A valid ISO3 country code.
                    
    """
    # modify the country attribute
    model.iso3_code = iso3_code
        
    # get the PWT data for the iso_code
    model.data      = model.pwt_data.minor_xs(iso3_code)
    model.dep_rates = model.pwt_dep_rates.minor_xs(iso3_code)
    
    ##### estimate capital's share of income/output ####
    alpha = (1 - model.data.labsh).mean()
            
    ##### estimate the fraction of output saved #####
    s = model.data.csh_i.mean()
        
    ##### estimate the labor force growth rate #####
        
    # regress log employed persons on linear time trend
    trend = pd.Series(model.data.index, index=model.data.index)
    res   = pd.ols(y=np.log(model.data.emp), x=trend)
    n     = res.beta[0]
    #L0    = np.exp(res.y_fitted[trend.index[0]])
    L0    = np.exp(res.beta[1])
    
    ##### estimate the technology growth rate #####
        
    # adjust measure of TFP to remove human capital contribution
    model.data['atfpna'] = model.data.rtfpna**(1 / (1 - alpha)) * model.data.hc
        
    # regress log TFP on linear time trend
    res   = pd.ols(y=np.log(model.data.atfpna), x=trend)
    g     = res.beta[0]
    #A0    = np.exp(res.y_fitted[trend.index[0]])
    A0    = np.exp(res.beta[1])
                           
    ##### estimate the depreciation rate for total capital #####
        
    # use average depreciation rate for total capital
    delta = model.dep_rates.delta_k.mean()   
         
    # create a dictionary of model parameters
    params = {'s':s, 'alpha':alpha, 'delta':delta, 'n':n, 'L0':L0, 'g':g,
              'A0':A0}
        
    # update the model's parameters
    model.update_model_parameters(params)
                    
    # compute new steady state values
    model.steady_state.set_values()

def calibrate_ces(model, iso3_code, x0, method='Nelder-Mead', tol=1e-9, 
                  bounds=None):
    """
    Calibrates a Solow model with constant elasticity of substition (CES)
    production using data from the Penn World Tables (PWT).

    Arguments:
        
        model:     (object) An instance of the SolowModel class.
            
        iso3_code: (str) A valid ISO3 country code.        
                   
        bounds:    (tuple) Start and end years for the subset of the PWT data
                   to use for calibration.
                   
    TODO: Validate the non-linear least squares approach being used to find
          optimal alpha and sigma.
        
    """
    # modify the country attribute
    model.iso3_code = iso3_code
        
    # get the PWT data for the iso_code
    model.data      = model.pwt_data.minor_xs(iso3_code)
    model.dep_rates = model.pwt_dep_rates.minor_xs(iso3_code)
    
    # set bounds
    if bounds == None:
        start = model.data.index[0]
        end   = model.data.index[-1]
    else:
        start = bounds[0]
        end   = bounds[1]
         
    ##### estimate the fraction of output saved #####
    s = model.data.csh_i.loc[start:end].mean()
        
    ##### estimate the labor force growth rate #####
        
    # regress log employed persons on linear time trend
    trend = pd.Series(model.data.index, index=model.data.index)
    res   = pd.ols(y=np.log(model.data.emp.loc[start:end]), 
                   x=trend.loc[start:end])
    n     = res.beta[0]
    #L0    = np.exp(res.y_fitted[trend.index[0]])
    L0    = np.exp(res.beta[1])
                           
    ##### estimate the depreciation rate for total capital #####
        
    # use average depreciation rate for total capital
    delta = model.dep_rates.delta_k.loc[start+1:end].mean()  
       
    ##### estimate alpha and sigma using NLLS #####
    
    def func(params):
        """Optimize to find alpha and sigma."""
        alpha = params[0]
        sigma = params[1]
        rho   = (sigma - 1) / sigma
    
        A = model.data.rtfpna.loc[start:end]
        Y = model.data.rgdpna.loc[start:end]
        K = model.data.rkna.loc[start:end]
        L = model.data.emp.loc[start:end]
        H = model.data.hc.loc[start:end]
    
        # adjust PWT TFP to A for CES production
        A_tilde = ((K / L) * (((A**rho * (K / L)**(-rho * (1 - alpha)) * 
                   H**(rho * (1 - alpha))) - alpha) / (1 - alpha))**(1 / rho)) 
    
        y = (Y / (A_tilde * L))
        k = (K / (A_tilde * L))
    
        # remove the discontinuity 
        if abs(rho) < 1e-6:
            return np.sum((np.log(y) - alpha * alpha * np.log(k))**2)
        else:
            return np.sum((np.log(y) - (1 / rho) * np.log(alpha * k**rho + (1 - alpha)))**2)
        
    # solve the non-linear least squares problem        
    res = optimize.minimize(func, x0, method=method, tol=tol)
    alpha, sigma = res.x
                
    ##### estimate the technology growth rate #####
          
    # adjust measure of TFP to remove human capital contribution
    A   = model.data.rtfpna
    K   = model.data.rkna
    L   = model.data.emp
    H   = model.data.hc
    rho = (sigma - 1) / sigma
    
    model.data['atfpna'] = ((K / L) * (((A**rho * (K / L)**(-rho * (1 - alpha)) * 
                            H**(rho * (1 - alpha))) - alpha) / (1 - alpha))**(1 / rho)) 

    # regress log TFP on linear time trend
    res   = pd.ols(y=np.log(model.data.atfpna.loc[start:end]), 
                   x=trend.loc[start:end])
                   
    g     = res.beta[0]
    #A0    = np.exp(res.y_fitted[trend.index[0]])
    A0    = np.exp(res.beta[1])
    
    # create a dictionary of model parameters
    params = {'s':s, 'alpha':alpha, 'sigma':sigma, 'delta':delta, 'n':n, 
              'L0':L0, 'g':g, 'A0':A0}
              
    # update the model's parameters
    model.update_model_parameters(params)
                    
    # compute new steady state values
    model.steady_state.set_values()
          
class RamseyModel(solvers.IVP):
    """Base class for a Ramsey (1928) model of optimal savings."""
    # each instance should carry a copy of the PWT data 
    pwt_data = pwt.load_pwt_data()
        
    def __init__(self, output, k_dot, c_dot, F, jacobian, params):
        """
        Initializes a RamseyModel object with the following attributes:
            
            output:   (callable) Output (per person/effective person). Should be
                      of the form
                      
                          f(t, k, params)
                          
                      where the independent variable, t, is time; k is capital
                      (per person/effective person); and params is a dictionary 
                      of model parameters.
                                             
            k_dot:    (callable) Equation of motion for capital (per person/
                      effective person). Should be of the form 
                     
                         k_dot(t, vec, params)
                        
                      where the independent variable, t, is time; vec is a vector
                      of the endogenous variables with ordering [k, c]; and 
                      params is a dictionary of model parameters.
                      
            c_dot:    (callable) Equation of motion for consumption (per person/
                      effective person). Should be of the form 
                     
                         c_dot(t, vec, params)
                        
                      where the independent variable, t, is time; vec is a vector
                      of the endogenous variables with ordering [k, c]; and 
                      params is a dictionary of model parameters.
                      
            F:        (callable) Function defining the 2D system of ODEs 
                      describing the evolution of the Ramsey economy. Should be
                      of the form
                      
                          F(t, vec, params)
                       
                      where the independent variable, t, is time; vec is a vector
                      of the endogenous variables with ordering [k, c]; and 
                      params is a dictionary of model parameters. 
            
            jacobian: (callable) Returns the Jacobian matrix of partial 
                      derivatives for F. Should be of the form
                      
                          jacobian(t, vec, params)
                       
                      where the independent variable, t, is time; vec is a vector
                      of the endogenous variables with ordering [k, c]; and 
                      params is a dictionary of model parameters.  
                                     
            params:   (dict) Dictionary of model parameters.
            
        """
        # initialize model attributes
        self.output                   = output
        self.k_dot                    = k_dot
        self.c_dot                    = c_dot
        self.iso3_code                = None
        self.data                     = None
        
        # initialize an empty SteadyState object
        self.steady_state = steady_states.SteadyState(self) 
        
        # initialize the model as an IVP
        super(RamseyModel, self).__init__(F, jacobian, params)            
            
    def __get_k_locus(self, t, k_grid, params):
        """Values of c consistent with steady state k."""
        # for each k, want value of c that makes this zero
        k_locus = lambda c, k: self.k_dot(0, [k, c], self.args)
        
        # Newton's method not guaranteed to converge!
        c_star = self.steady_state.values['c_star']
        out = [optimize.newton(k_locus, c_star, args=(k,)) for k in k_grid]
        
        return np.array(out)
        
    def plot_phase_diagram(self, gridmax, N=1000, arrows=False, param=None,  
                           shock=None, reset=True, cmap='winter', mu=0.1):
        """
        Generates phase diagram for the Ramsey model.
        
        Arguments:
            
            gridmax: (float) Maximum value for capital per effective person.
            
            N:       (int, optional) Number of points to plot. Default is 1000.
            
            arrows:  (boolean) If True, plots directional arrows indicating out
                     of steady state dynamics. Default is False.
            
            param:   (string) Model parameter to shock (optional).
            
            shock:   (float) Multiplicative shock to the parameter. Values < 1 
                     correspond to a reducing the current value of a parameter; 
                     values > 1 correspond to increasing the current value of 
                     the parameter (optional).
            
            reset:   (boolean) Whether or not to reset the original parameters
                     to their pre-shock values. Default is True.
            
            cmap:    (str) A valid matplotlib colormap. Default is 'winter'.
            
            mu:      (float) Determines spread between colors...
        
        Returns: A list containing...
        
            ax:          (object) Axes object representing the plot.
            
            k_locus:     (object) Line2D object representing the k-dot locus.
            
            c_locus:     (object) Line2D object representing actual c-dot locus.
            
            ss_marker:   (object) Line2D object representing the steady state.
            
        """
        # sets the color palette
        colors = mpl.cm.__dict__[cmap]([mu, (1 - mu)])
               
        # create a new figure and subplot
        ax     = plt.subplot(111)

        # use steady state values to anchor the plot
        k_star = self.steady_state.values['k_star']
        c_star = self.steady_state.values['c_star']
        
        # grid of points for plotting
        grid   = np.linspace(0, gridmax, N)
        
        k_locus   = ax.plot(grid, self.__get_k_locus(0, grid, self.args), '-', 
                            color=colors[0], label=r'$\dot{k}=0$')[0]
        
        c_locus   = ax.axvline(k_star, color=colors[1], label=r'$\dot{c}=0$')
        
        ss_marker = ax.plot(k_star, c_star, marker='.', markersize=10, 
                            color='k')[0]

        # remove the right and top spines
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')

        # hide the top and right ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # demarcate the steady state
        ax.set_xticks([k_star])
        ax.set_xticklabels([r'$k^*$'])
        ax.set_yticks([c_star])
        ax.set_yticklabels([r'$c^*$'])
        
        # axes, labels, title, legend, etc
        ax.set_xlabel('$k_t$', fontsize=15)
        ax.set_ylim(0, 2 * c_star)
        ax.set_ylabel('$c_t$', rotation='horizontal', fontsize=15)
        
        # Add arrows to indicate out of steady-state dynamics
        if arrows == True:
            
            mu = 0.25
            x_len = mu * k_star 
            y_len = mu * c_star   

            ax.arrow(x=0.5 * k_star, y=0.5 * c_star, dx=0, dy=y_len, fc='k',
                     shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * k_star, head_length=0.05 * c_star)
            ax.arrow(x=0.5 * k_star, y=0.5 * c_star, dx=x_len, dy=0, fc='k',
                     shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * c_star, head_length=0.05 * k_star)

            ax.arrow(x=0.5 * k_star + x_len, y=1.5 * c_star, dx=0, dy=y_len,
                     fc='k', shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * k_star, head_length=0.05 * c_star)
            ax.arrow(x=0.5 * k_star + x_len, y=1.5 * c_star, dx=-x_len, dy=0,
                     fc='k', shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * c_star, head_length=0.05 * k_star)

            ax.arrow(x=1.5 * k_star, y=0.5 * c_star + y_len, dx=0, dy=-y_len,
                     fc='k', shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * k_star, head_length=0.05 * c_star)
            ax.arrow(x=1.5 * k_star, y=0.5 * c_star + y_len, dx=x_len, dy=0,
                     fc='k', shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * c_star, head_length=0.05 * k_star)

            ax.arrow(x=1.5 * k_star + x_len, y=1.5 * c_star + y_len, dx=0, 
                     dy=-y_len, fc='k', shape='full', lw=1, 
                     length_includes_head=True, head_width=0.05 * k_star, 
                     head_length=0.05 * c_star)
            ax.arrow(x=1.5 * k_star + x_len, y=1.5 * c_star + y_len, dx=-x_len, 
                     dy=0, fc='k', shape='full', lw=1, length_includes_head=True, 
                     head_width=0.05 * c_star, head_length=0.05 * k_star)

        # handles parameter shocks   
        if param != None and shock !=None:
            
            # copy the original params
            orig_params = self.args.copy()
            
            # shock the parameter
            self.args[param] = shock * self.args[param]
        
            # compute the new steady state values
            self.steady_state.set_values()
            k_star = self.steady_state.values['k_star']
            c_star = self.steady_state.values['c_star']
        
            # demarcate the new steady state
            new_ss_marker = ax.plot(k_star, c_star, marker='.', markersize=10, 
                                    color='k')[0]
            ax.set_xticks([k_star])
            ax.set_xticklabels([r'$k^*$'])
            ax.set_yticks([c_star])
            ax.set_yticklabels([r'$c^*$'])
            
            # reset y-axes limits
            ax.set_ylim(0, 2 * c_star)
            
            # plot formatting depends on parameter being shocked
            if param in ['alpha', 'delta', 'rho', 'theta']:
                param = '\\' + param  # necessary for pretty latex printing
                
            # changes in these params shift both k and c locii
            if param in ['\\alpha', 'g', '\\delta']:
                new_k_locus = ax.plot(grid, self.__get_k_locus(0, grid, self.args), 
                                      color=colors[0], label=r'$\dot{k}=0$')[0]
                new_c_locus = ax.axvline(k_star, color=colors[1], 
                                         label=r'$\dot{c}=0$')
                
                k_locus.set_alpha(0.5)
                k_locus.set_linestyle('dashed')
                k_locus.set_label('$\dot{k}=0_{old}$')
                new_k_locus.set_label('$\dot{k}=0_{new}$')
                
                c_locus.set_alpha(0.5)
                c_locus.set_linestyle('dashed')
                c_locus.set_label('$\dot{c}=0_{old}$')
                new_c_locus.set_label('$\dot{c}=0_{new}$')
                
                ss_marker.set_alpha(0.5)
                
                ax.set_title(('Changing $%s$ shifts both $\dot{k}=0$ and ' + 
                              '$\dot{c}=0$ locii!') %param, fontsize=20, 
                              family='serif')
                ax.legend(loc='best', frameon=False)
                
            # changes in these params shift the c-dot locus only
            elif param in ['\\rho', '\\theta']:
                new_k_locus = None
                new_c_locus = ax.axvline(k_star, color=colors[1], 
                                         label=r'$\dot{c}=0$')
                c_locus.set_alpha(0.5)
                c_locus.set_linestyle('dashed')
                c_locus.set_label('$\dot{c}=0_{old}$')
                new_c_locus.set_label('$\dot{c}=0_{new}$')
                
                ss_marker.set_alpha(0.5)
                
                ax.set_title('Changing $%s$ shifts the $\dot{c}=0$ locus!' 
                              %param, fontsize=20, family='serif')
                ax.legend(loc='best', frameon=False) 
            
            # changes in the population growth rate shift k-dot locus
            elif param == 'n':
                new_k_locus = ax.plot(grid, self.__get_k_locus(0, grid, self.args), 
                                      color=colors[0], label=r'$\dot{k}=0$')[0]
                new_c_locus = None
                
                k_locus.set_alpha(0.5)
                k_locus.set_linestyle('dashed')
                k_locus.set_label('$\dot{k}=0_{old}$')
                new_k_locus.set_label('$\dot{k}=0_{new}$')
                
                ss_marker.set_alpha(0.5)
                
                ax.set_title('Changing $%s$ shifts the $\dot{k}=0$ locus!' 
                              %param, fontsize=20, family='serif')
                ax.legend(loc='best', frameon=False)
                
            else:
                raise ValueError
            
            # reset the original params and recompute steady state?
            if reset == True:
                self.args = orig_params
                self.steady_state.set_values()
                
            return [ax, new_k_locus, new_c_locus, new_ss_marker]
        
        else:
            ax.set_title('Phase diagram for the Ramsey (1928) model', fontsize=20, 
                         family='serif')
            ax.legend(loc='best', frameon=False)
        
            return [ax, k_locus, c_locus, ss_marker]   
        
    def solve_forward_shooting(self, k0, h=1e-3, tol=1.5e-3, mesg=False, 
                               integrator='lsoda', **kwargs):
        """
        Computes the full, non-linear saddle path for the continuous time 
        version of the Ramsey model using the 'forward shooting' algorithm (see 
        Judd (1998) p. 357 for details).

        Arguments:

            k0:         (float) Initial value for capital (per person/effective
                        person).
            
            h:          (float) Step-size for computing the solution trajectory.
            
            tol:        (float) Convergence tolerance for solution trajectory. 
                        Due to accumulation of numerical error in solving the 
                        IVP for some k0 and c0, it may be necessary to choose a
                        relatively loose stopping criterion.
                 
            mesg:       (boolean) If True, then messages are printed showing 
                        convergence progress. Default is False.
                         
            integrator: (str) Must be a valid integrator class.See docstring of 
                        the integrate method for a complete listing of valid
                        integrators. 
                     
            **kwargs:   (dict) Dictionary of integrator specific keyword args.
                
        Returns: 
                     
            solution: (array-like) Simulated solution trajectory approximating 
                      the model's saddle-path.
               
        """ 
        # compute steady state values
        k_star = self.steady_state.values['k_star']
        c_star = self.steady_state.values['c_star']
   
        # compute the bounds for initial guess
        if k0 < k_star:
            c_l = 0
            c_h = self.__get_k_locus(0, [k0], self.args)[0]
        else:
            c_l = self.__get_k_locus(0, [k0], self.args)[0]
            c_h = (1 - self.args['delta']) * k0 + self.output(0, k0, self.args)

        # default initial guess for c 
        c0 = (c_h + c_l) / 2
        
        # create an instance of the scipy.integrate.ode class      
        ode  = integrate.ode(self.f, self.jac)
        
        # select the integrator
        ode.set_integrator(integrator, **kwargs)
        
        # pass the model parameters as additional args
        ode.set_f_params(self.args)
        ode.set_jac_params(self.args)
        
        # set the initial condition
        ode.set_initial_value([k0, c0], 0)
        
        ########## Compute the optimal c0 using bisection method ############
        while ode.successful():
            # integrate the system one step
            ode.integrate(ode.t + h)
            
            # get the values of the vector field
            k_dot, c_dot = self.f(ode.t, ode.y, self.args)
            
            if k0 < k_star and k_dot < 0:
                
                # check for convergence
                if abs(ode.y[1] - c_star) < tol:
                    ode.set_initial_value([k0, c0], 0)
                    break
                
                # c0 too high!
                else:
                    c_h = c0
                    c0  = 0.5 * (c_h + c_l)
                    if mesg == True:
                        print 'Old c0 too high, new c0 =', c0
                    ode.set_initial_value([k0, c0], 0)
                      
            elif k0 < k_star and c_dot < 0:    
                
                # check for convergence
                if abs(ode.y[1] - c_star) < tol:
                    ode.set_initial_value([k0, c0], 0)
                    break
                    
                # c0 too low!
                else:
                    c_l = c0
                    c0  = 0.5 * (c_h + c_l)
                    if mesg == True:
                        print 'Old c0 too low, new c0 =', c0
                    ode.set_initial_value([k0, c0], 0)
                    
            elif k0 > k_star and k_dot > 0:
                
                # check for convergence
                if abs(ode.y[1] - c_star) < tol:
                    ode.set_initial_value([k0, c0], 0)
                    break
                
                # c0 too low!
                else:
                    c_l = c0
                    c0  = 0.5 * (c_h + c_l)
                    if mesg == True:
                        print 'Old c0 too low, new c0 =', c0
                    ode.set_initial_value([k0, c0], 0)
                      
            elif k0 > k_star and c_dot > 0:    
                
                # check for convergence
                if abs(ode.y[1] - c_star) < tol:
                    ode.set_initial_value([k0, c0], 0)
                    break
                    
                # c0 too high!
                else:
                    c_h = c0
                    c0  = 0.5 * (c_h + c_l)
                    if mesg == True:
                        print 'Old c0 too high, new c0 =', c0
                    ode.set_initial_value([k0, c0], 0)
            
            else:
                continue
        
        ########## Compute the saddle path using the optimal c0 ##########
                
        # create a storage container for the trajectory
        solution = np.hstack((0, ode.y)) 
          
        while ode.successful() and abs(ode.y[1] - c_star) > tol:
            ode.integrate(ode.t + h)
            
            # store the current step
            current_step = np.hstack((ode.t, ode.y))
            solution = np.vstack((solution, current_step))  
                
        return solution
        
    def solve_reverse_shooting(self, k0, h=1e-3, eps=1e-5, integrator='dopri5', 
                               step=False, relax=False, **kwargs):
        """
        Computes the full, non-linear saddle path (i.e., the consumption policy
        function) using a 'reverse shooting' algorithm (see Judd (1992) section 
        10.7 Infinite-Horizon Optimal Control and Reverse Shooting, p. 355-361 
        for details).
        
        Arguments:
                            
            k0:         (float) Initial condition for capital (per person/
                        effective person)
                                                                 
            h:          (float) Step-size for computing the solution trajectory.
            
            eps:        (float) Initial step size.
            
            integrator: (str) Must be a valid integrator class.See docstring of 
                        the integrate method for a complete listing of valid
                        integrators. Default is 'dopri5'.
                     
            **kwargs:   (dict) Dictionary of integrator specific keyword args.
                
        Returns: 
                     
           solution: (array-like) Simulated solution trajectory approximating 
                      the model's saddle-path.
               
        """ 
        # compute steady state values
        k_star = self.steady_state.values['k_star']
        c_star = self.steady_state.values['c_star']

        # find index of the stable eigenvalue
        index = np.where(np.real(self.steady_state.eigenvalues) < 0)[0][0]
        
        # local slope of optimal policy evaluated at steady state
        Ck_prime = (self.steady_state.eigenvectors[1, index] / 
                    self.steady_state.eigenvectors[0, index])
        
        # RHS of equation 10.7.5 from Judd (1998)
        c_prime = lambda k, c, params: (self.c_dot(0, [k, c], params) / 
                                        self.k_dot(0, [k, c], params))
         
        # create an instance of the scipy.integrate.ode class      
        ode = integrate.ode(c_prime)
        
        # select the integrator
        ode.set_integrator(integrator, **kwargs)
        
        # pass the model parameters as additional args
        ode.set_f_params(self.args)
        
        # set initial conditions
        if k0 > k_star:
            ode.set_initial_value(c_star + eps * Ck_prime, k_star + eps)
        else:
            ode.set_initial_value(c_star - eps * Ck_prime, k_star - eps)
        
        # create a storage container for the trajectory
        solution = np.hstack((ode.t, ode.y)) 
        
        # generate a solution trajectory
        while ode.successful():               
            ode.integrate(ode.t + h, step, relax)
            current_step = np.hstack((ode.t, ode.y))
            solution = np.vstack((solution, current_step))  
            
            # check to see if initial condition has been reached
            if k0 < k_star and ode.t < k0:
                break
            elif k0 > k_star and ode.t > k0:
                break
            else:
                continue 
        
        return solution
    
    def solve_multiple_shooting(self, k0, T, solution_guess, 
                                boundary_conditions, **kwargs):
        """
        Wraps scikits.bvp_solver in order to solve the model using a multiple
        shooting approach.
        
        Arguments:
            
            k0:                  (float) Initial condition for capital (per 
                                 person/effective person)
                        
            T:                   (float) Right boundary of the time interval of 
                                 interest Needs to be large enough to ensure  
                                 that the system will have enough time to 
                                 converge to steady state.
                            
            solution_guess:      (Solution, float, array) An initial guess for
                                 the true solution. 
            
            boundary_conditions: (callable): A function with calculates the 
                                 difference between the actual boundary 
                                 conditions and the desired boundary conditions.
                
            **kwargs:            (dict) Dictionary of keyword arguments passed 
                                 to the scikits.bvp_solver.solve method.
                            
        Returns:
            
            result: (object) An instance of the scikits.bvp_solver.Solution 
                    class.
            
        """
        # scikits.bvp_solver requires slightly different f and jac
        f   = lambda t, vec: self.f(t, vec, self.args)
        jac = lambda t, vec: self.jac(t, vec, self.args)
                
        # Create the ProblemDefinition argument
        bvp = bvp_solver.ProblemDefinition(num_ODE = 2,
                                           num_parameters = 0,
                                           num_left_boundary_conditions = 1,
                                           boundary_points = (0, T),
                                           function = f,
                                           function_derivative = jac,
                                           boundary_conditions = boundary_conditions)
        
        sol = bvp_solver.solve(bvp, solution_guess = solution_guess, **kwargs)
        
        return sol
        
    def get_stable_manifold(self, kmin, kmax, method=None, **kwargs):
        """
        Computes the stable manifold for the Ramsey model using either 'forward'
        or 'reverse' shooting, depending on whether 'tol' or 'eps' is specified.
        
        Arguments:
                        
            kmin:       (float) Terminal condition for capital (per person/
                        effective person) for the lower portion of the stable
                        manifold.
                                    
            kmax:   	(float) Terminal condition for capital (per person/
                        effective person) for the upper portion of the stable
                        manifold.
            
            method:     (str) One of 'forward', 'reverse', or 'multiple', 
                        depending.
                        
            **kwargs:   (dict) Dictionary of method specific keyword args. For
                        method = 'forward' the following keyword arguments are 
                        required:
                            
                            h:          (float) Step-size to use for the 
                                        integration. Default is 1.0.
                                        
                            tol:        (float) Convergence tolerance for 
                                        solution trajectory. Default is 1e-6.
                                        
                            integrator: (str) Must be a valid integrator class. 
                                        See docstring of the integrate method 
                                        for a complete listing of valid
                                        integrators. Default is 'dopri'.
                                        
                            options:    (dict) Dictionary of integrator specific
                                        keyword arguments. Default is {}.
                                                                               
                        For method = 'reverse' the following keyword arguments 
                        are required:
                            
                            h:          (float) Step-size to use for the 
                                        integration. Default is 1.0.
                                        
                            eps:        (float) Initial step-size. Default is 
                                        1e-6.
                                        
                            integrator: (str) Must be a valid integrator class. 
                                        See docstring of the integrate method 
                                        for a complete listing of valid
                                        integrators. Default is 'dopri'.
                                        
                            options:    (dict) Dictionary of integrator specific
                                        keyword arguments. Default is {}.
                                        
                        For method = 'multiple' the following keyword arguments
                        are required:
                            
                            T:                   (float) Upper boundary for the 
                                                 time interval of interest.
                                                 Default is 1000.
                                                 
                            solution_guess:      (Solution, float, array) An 
                                                 initial guess for the true 
                                                 solution. Default is None.
            
                            boundary_conditions: (callable): A function with 
                                                 calculates the difference 
                                                 between the actual boundary 
                                                 conditions and the desired
                                                 boundary conditions. Default is
                                                 None.
                            
                            options:             (dict) Dictionary of keyword 
                                                 arguments to pass to the 
                                                 scikits.bvp_solver.solve 
                                                 method. Default is {}.

        Returns:
            
            M_S: (array-like) Array representing the stable manifold for the 
                 Ramsey model. 
                 
        """
        if method == 'reverse':
            
            # extract the keyword args for 'reverse' shooting
            h           = kwargs.get('h', 1.0)
            eps         = kwargs.get('eps', 1e-6)
            integrator  = kwargs.get('integrator', 'dopri')
            options     = kwargs.get('options', {})
            
            lower_M_S = self.solve_reverse_shooting(kmin, -h, eps, integrator, 
                                                    **options)
            upper_M_S = self.solve_reverse_shooting(kmax, h, eps, integrator, 
                                                    **options)
          
            # reverse the direction of lower_M_S
            lower_M_S = lower_M_S[::-1]
        
        elif method == 'forward':
            
            # extract the keyword args for 'forward' shooting
            h           = kwargs.get('h', 1.0)
            tol         = kwargs.get('tol', 1e-6)
            integrator  = kwargs.get('integrator', 'dopri')
            options     = kwargs.get('options', {})
            
            lower_M_S = self.solve_forward_shooting(kmin, h, tol, integrator, 
                                                    **options)
            upper_M_S = self.solve_forward_shooting(kmax, h, tol, integrator, 
                                                    **options)
            
            # reverse the direction of upper_M_S
            upper_M_S = upper_M_S[::-1]
            
            # drop the time dimension
            lower_M_S = lower_M_S[:,1:]
            upper_M_S = upper_M_S[:,1:]
            
        elif method == 'multiple':
            
            # extract the keywords for 'multiple' shooting
            T                   = kwargs.get('T', 1e3)
            solution_guess      = kwargs.get('solution_guess', None)
            boundary_conditions = kwargs.get('boundary_conditions', None)
            options             = kwargs.get('options', {})
            c_star              = self.steady_state.values['c_star']
            
            # wrap the boundary conditions for k0 = kmin
            bc = lambda veca, vecb: boundary_conditions(veca, vecb, kmin, c_star)
            lower = self.solve_multiple_shooting(kmin, T, solution_guess, bc, 
                                                 **options)
                       
            # wrap the boundary conditions for k0 = kmax                          
            bc = lambda veca, vecb: boundary_conditions(veca, vecb, kmax, c_star)
            upper = self.solve_multiple_shooting(kmax, T, solution_guess, 
                                                 bc, **options)
                     
            # bvp_solver returns trajectories as row vectors!
            lower_M_S = lower.solution.T
            upper_M_S = upper.solution.T
            
            # reverse the direction of upper_M_S
            upper_M_S = upper_M_S[::-1]
                              
        else:
            raise ValueError 
            
        # average the last row of lower_M_S with first of upper_M_S
        lower_M_S[-1] = 0.5 * (lower_M_S[-1] + upper_M_S[0])
        
        # combine the lower and upper portions of M_S
        M_S = np.vstack((lower_M_S, upper_M_S[1:]))
        
        return M_S

    def get_consumption_policy(self, kmin, kmax, method='multiple', **kwargs):
        """
        Constructs a callable representation of the representative household's
        consumption policy function. 
        
        Arguments:
                        
            kmin:       (float) Terminal condition for capital (per person/
                        effective person) for the lower portion of the stable
                        manifold.
                                    
            kmax:   	(float) Terminal condition for capital (per person/
                        effective person) for the upper portion of the stable
                        manifold.
                                 
            method:     (str) One of 'forward', 'reverse', or 'multiple', 
                        depending.
                        
            **kwargs:   (dict) Dictionary of method specific keyword args. For
                        method = 'forward' the following keyword arguments are 
                        required:
                            
                            h:          (float) Step-size to use for the 
                                        integration. Default is 1.0.
                                        
                            tol:        (float) Convergence tolerance for 
                                        solution trajectory. Default is 1e-6.
                                        
                            integrator: (str) Must be a valid integrator class. 
                                        See docstring of the integrate method 
                                        for a complete listing of valid
                                        integrators. Default is 'dopri'.
                                        
                            options:    (dict) Dictionary of integrator specific
                                        keyword arguments. Default is {}.
                                            
                            k:          (int) Degree of the desired B-spline. 
                                        Must satisfy 1 <= k <= 5. Default is 3.
                                        
                        For method = 'reverse' the following keyword arguments 
                        are required:
                            
                            h:          (float) Step-size to use for the 
                                        integration. Default is 1.0.
                                        
                            eps:        (float) Initial step size. Default is 
                                        1e-6.
                                        
                            integrator: (str) Must be a valid integrator class. 
                                        See docstring of the integrate method 
                                        for a complete listing of valid
                                        integrators. Default is 'dopri'.
                                        
                            options:    (dict) Dictionary of integrator specific
                                        keyword arguments. Default is {}.
                                        
                            k:          (int) Degree of the desired B-spline. 
                                        Must satisfy 1 <= k <= 5. Default is 3.
                                        
                       For method = 'multiple' the following keyword arguments
                       are required:
                            
                            T:                   (float) Upper boundary for the 
                                                 time interval of interest.
                                                 Default is 1000.
                                                 
                            solution_guess:      (Solution, float, array) An 
                                                 initial guess for the true 
                                                 solution. Default is None.
            
                            boundary_conditions: (callable): A function with 
                                                 calculates the difference 
                                                 between the actual boundary 
                                                 conditions and the desired
                                                 boundary conditions. Default is
                                                 None.
                            
                            options:             (dict) Dictionary of keyword 
                                                 arguments to pass to the 
                                                 scikits.bvp_solver.solve 
                                                 method. Default is {}.

        Returns:
            
            c_k: (callable) A callable object representing the consumption 
                 policy function.
                 
        """
        # extract keyword args
        k = kwargs.get('k', 3)
            
        # compute the stable manifold
        M_S = self.get_stable_manifold(kmin, kmax, method, **kwargs)
        
        # construct a callable B-spline repr of the policy function
        c_k = interpolate.UnivariateSpline(M_S[:,0], M_S[:,1], k=k, s=0)                   
                            
        return c_k

    def compare_policies(self, pol1, pol2, grid):
        """
        Compare two consumption policy functions at common grid points.
        
        Arguments:
            
            pol1: (callable) Consumption policy function computed using the 
                  get_consumption_policy method.
                  
            pol2: (callable) Consumption policy function computed using the 
                  get_consumption_policy method.
                  
            grid: (array-like) Grid of values for capital per effective worker
                  at which to compare the two policy functions.
                  
        Returns:
            
            diff: (array) Array of element-wise differences between the two 
                  policy functions.
        
        """
        diff = pol1(grid) - pol2(grid)
        return diff