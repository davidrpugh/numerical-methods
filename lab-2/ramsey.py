from __future__ import division

import pandas as pd
import numpy as np
from scipy import integrate, interpolate, linalg, optimize
#import scikits.bvp_solver as bvp_solver

import matplotlib as mpl
import matplotlib.pyplot as plt

import solvers, integrators
import steady_states
import pwt
                    
class Model(solvers.IVP):
    """Base class for a Ramsey (1928) model of optimal savings."""
    # each instance should carry a copy of the PWT data 
    pwt_data, pwt_dep_rates = pwt.load_pwt_data(deltas=True)
        
    def __init__(self, output, mpk, k_dot, c_dot, utility, jacobian, params):
        """
        Initializes a RamseyModel object with the following attributes:
            
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
                     
                         k_dot(t, vec, params)
                        
                      where the independent variable, t, is time; vec is a
                      vector of the endogenous variables with ordering [k, c]; 
                      and params is a dictionary of model parameters.
                      
            c_dot:    (callable) Equation of motion for consumption (per person/
                      effective person). Should be of the form 
                     
                         c_dot(t, vec, params)
                        
                      where the independent variable, t, is time; vec is a 
                      vector of the endogenous variables with ordering [k, c]; 
                      and params is a dictionary of model parameters.
                      
            utility:  (callable) Function defining the instantaneous utility 
                      function used to derive the consumption Euler equation.
                      Should be of the form 
                     
                         u(c, params)
                        
                      where c is consumption (per person/effective person) and 
                      params is a dictionary of model parameters.
            
            jacobian: (callable) Returns the Jacobian matrix of partial 
                      derivatives for the model. Should be of the form
                      
                          jacobian(t, vec, params)
                       
                      where the independent variable, t, is time; vec is a
                      vector of the endogenous variables with ordering [k, c]; 
                      and params is a dictionary of model parameters.  
                                     
            params:   (dict) Dictionary of model parameters.
            
        """
        # initialize model attributes
        self.output                   = output
        self.mpk                      = mpk
        self.k_dot                    = k_dot
        self.c_dot                    = c_dot
        self.iso3_code                = None
        self.data                     = None
        self.dep_rates                = None
        
        # combine k_dot and c_dot into a system of equations
        def F(t, vec, params):
            """System of equations defining the optimal growth model."""
            out = [self.k_dot(t, vec, params), self.c_dot(t, vec, params)]
            return np.array(out)
           
        # initialize the model as an IVP 
        super(Model, self).__init__(F, jacobian, (params,))
        self.params = self.args[0]
        
        # initialize an empty SteadyState object
        self.steady_state = steady_states.SteadyState(self)    
                            
    def update_model_parameters(self, new_params):
        """Updates the model's parameter dictionary."""
        self.args = (new_params.copy(),) + self.args[1:]
        self.params = self.args[0]            
    
    def get_capitals_share(self, k):
        """Computes capital's share of income eps_k(k)."""
        capitals_income = k * self.mpk(0, k, self.params)
        total_output    = self.output(0, k, self.params)
        
        # definition of capital's share
        eps_k =  capitals_income / total_output 
    
        return eps_k
                
    def __get_k_locus(self, t, k_grid, params):
        """Values of c consistent with steady state k."""
        # for each k, want value of c that makes this zero
        k_locus = lambda c, k: self.k_dot(0, [k, c], self.params)
        
        # Newton's method not necessarily guaranteed to converge!
        c_star = self.steady_state.values['c_star']
        out    = [optimize.newton(k_locus, c_star, args=(k,)) for k in k_grid]
        
        return np.array(out)
        
    def __add_arrows(self, ax, mu=0.25):
        """
        Modifies the phase diagram by adding directional arrows.
        
        Arguments:
            
            ax: (object) AxesSubplot object containing a phase diagram.
            
            mu: (float) Controls size of arrows.
            
        """
        # use steady state values to anchor the plot
        k_star = self.steady_state.values['k_star']
        c_star = self.steady_state.values['c_star']
        
        # define arrow size
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
                     
    def plot_phase_diagram(self, gridmax, N=1000, arrows=False, param=None,  
                           shock=None, reset=True, plot_traj=False, 
                           cmap='winter', mu=0.1):
        """
        Generates phase diagram for the Ramsey model.
        
        Arguments:
            
            gridmax:   (float) Maximum value for capital per effective person.
            
            N:         (int) Number of points to plot. Default is 1000.
            
            arrows:    (boolean) If True, plots directional arrows indicating
                       out of steady state dynamics. Default is False.
            
            param:     (string) Model parameter to shock (optional).
            
            shock:     (float) Multiplicative shock to the parameter. Values < 1 
                       correspond to a reducing the current value of a parameter; 
                       values > 1 correspond to increasing the current value of 
                       the parameter (optional).
            
            reset:     (boolean) Whether or not to reset the original parameters
                       to their pre-shock values. Default is True.
            
            plot_traj: (boolean) Whether or not you wish to plot a trajectory of
                       the economy in order to show transition dynamics 
                       following a shock to a parameter.
                       
            cmap:      (str) A valid matplotlib colormap. Default is 'winter'.
            
            mu:        (float) Determines spread between colors...
        
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
        
        # demarcate model steady state
        k_locus   = ax.plot(grid, self.__get_k_locus(0, grid, self.params), '-', 
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
        ax.set_xlabel('$k$', fontsize=15)
        ax.xaxis.set_label_coords(0.95, -0.05)
        ax.set_ylim(0, 2 * c_star)
        ax.set_ylabel('$c$', rotation='horizontal', fontsize=15)    
        ax.yaxis.set_label_coords(-0.05, 0.95)
        
        # handles parameter shocks   
        if param != None and shock !=None:
            
            # copy the original params and steady state value
            orig_params = self.params.copy()
            orig_k_star = k_star
            
            # shock the parameter
            self.params[param] = shock * self.params[param]
            self.update_model_parameters(self.params)
        
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
            if param in ['alpha', 'delta', 'rho', 'theta', 'sigma']:
                param = '\\' + param  # necessary for pretty latex printing
                
            # changes in these params shift both k and c locii
            if param in ['\\alpha', 'g', '\\delta', '\\sigma']:
                new_k_locus = ax.plot(grid, self.__get_k_locus(0, grid, self.params), 
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
                new_k_locus = ax.plot(grid, self.__get_k_locus(0, grid, self.params), 
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
            
            # plot a trajectory in phase space showing transition dynamics
            if plot_traj == True:
                traj = self.solve_reverse_shooting(orig_k_star)
                self.plot_trajectory(traj, color='r')
                
            # add arrows to indicate out of steady-state dynamics?
            if arrows == True:
                self.__add_arrows(ax, mu=0.25)
                
            # reset the original params and recompute steady state?
            if reset == True:
                self.update_model_parameters(orig_params)
                self.steady_state.set_values()
                
            return [ax, new_k_locus, new_c_locus, new_ss_marker]
        
        else:
            ax.set_title('Phase diagram for the optimal growth model', 
                         fontsize=20, family='serif')
            ax.legend(loc='best', frameon=False)
        
            # Add arrows to indicate out of steady-state dynamics?
            if arrows == True:
                self.__add_arrows(ax, mu=0.25)
                
            return [ax, k_locus, c_locus, ss_marker]   
        
    def plot_trajectory(self, traj, color='b', ax=None, phase_space=True, 
                        **kwargs):
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
                   
            kind:     (str) One of 'phase_space' or 'time_series', depending.
            
            **kwargs: (dict) Dictionary of optional keyword arguments to pass
                      to ax.plot method.
                   
        """ 
        if ax == None:
            ax = plt.subplot(111)
                
        if phase_space == True:
            ax.plot(traj[:,1], traj[:,2], color=color, **kwargs)
            
        else:
            ax.plot(traj[:,0], traj[:,1], color=color, **kwargs)
    
    def get_impulse_response(self, method, param, shock, T=100, 
                             kind='efficiency_units', reset=True):
        """
        Generates an impulse response functions for the endogenous variables in
        the model following a shock to one of the model parameters.
        
        Arguments:
            
            method: (str) One of 'linearization', or 'forward_shooting'.
            
            param: (string) Model parameter that you wish to shock.
            
            shock: (float) Multiplicative shock to the parameter. Values < 1 
                   correspond to a reduction in the current value of the 
                   parameter; values > 1 correspond to increasing the current 
                   value of the parameter.
                   
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
        c0 = self.steady_state.values['c_star']
        y0 = self.output(0, k0, self.params)
        r0 = self.mpk(0, k0, self.params) - self.params['delta']
        w0 = y0 - k0 * r0
        
        # initial padding should be such that shock occurs in '2013'
        N = 2013 - self.data.index[0] + 1
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
        
        # padding for c
        padding_c = np.repeat(c0, N)
        padding = np.hstack((padding, (factor * padding_c)[:,np.newaxis]))
        
        # padding for y
        padding_y = np.repeat(y0, N)
        padding = np.hstack((padding, (factor * padding_y)[:,np.newaxis]))
        
        # padding for r (r is same regardless of 'kind')
        padding_r = np.repeat(r0, N)
        padding = np.hstack((padding, padding_r[:,np.newaxis]))
        
        # padding for w 
        padding_w = np.repeat(w0, N)
        padding = np.hstack((padding, (factor * padding_w)[:,np.newaxis]))
        
        # move time padding forward to start of available data
        padding[:,0] += self.data.index[0]
                    
        # shock the parameter
        self.params[param] = shock * self.params[param]
        self.update_model_parameters(self.params)
        
        # compute the new steady state values
        self.steady_state.set_values()
        
        # generate post-shock trajectory 
        if method == 'linearization':
            ti  = np.linspace(0, T, T)
            irf = self.solve_linearization(k0, ti) 
                  
        elif method == 'forward_shooting':
            irf = self.solve_forward_shooting(k0)
            
            # compute some extra repeats of steady state values for padding
            k_star = self.steady_state.values['k_star']
            c_star = self.steady_state.values['c_star']
            
            extra_time = np.arange(int(irf[-1,0]), T, 1)[:,np.newaxis]
            extra_kss  = np.repeat(k_star, T - int(irf[-1,0]))[:,np.newaxis]
            extra_css  = np.repeat(c_star, T - int(irf[-1,0]))[:,np.newaxis]
            extra_padding = np.hstack((extra_time, extra_kss, extra_css))

            irf = np.vstack((irf, extra_padding))
            
        else:
            raise ValueError, "Invalid 'method' for computing irfs."    
        
        # transform irfs into per capita or levels, depending
        if kind == 'per_capita':
            g      = self.params['g']            
            factor = factor[-1] * np.exp(g * irf[:,0])
            
        elif kind == 'levels':
            g      = self.params['g']
            n      = self.params['n'] 
            factor = factor[-1] * np.exp((n + g) * irf[:,0])
            
        elif kind == 'efficiency_units':
            factor = np.ones(irf[:,0].size)
        
        else:
            raise ValueError, "Invalid 'kind' of impulse response functions."
        
        # shift the time index forward 2013 periods
        irf[:,0] += 2013
        
        # compute the irf for c 
        irf[:,2] = (factor * irf[:,2])
        
        # compute the irf for y
        irf_y = self.output(irf[:,0], irf[:,1], self.params)
        irf = np.hstack((irf, (factor * irf_y)[:,np.newaxis]))
        
        # compute the irf for r
        irf_r = self.mpk(irf[:,0], irf[:,1], self.params) - self.params['delta']
        irf = np.hstack((irf, irf_r[:,np.newaxis]))
        
        # compute the irf for w
        irf_w = irf_y - irf_r * irf[:,1]
        irf = np.hstack((irf, (factor * irf_w)[:,np.newaxis]))
        
        # compute the irf for k (after computing irfs for y, r, and w!)
        irf[:,1] = (factor * irf[:,1])
        
        # add the padding
        irf = np.vstack((padding, irf))
        
        # reset the original params and recompute steady state?
        if reset == True:
            self.update_model_parameters(orig_params)
            self.steady_state.set_values()
        
        return irf 
    
    def plot_impulse_response(self, variables, method, param, shock, T, 
                              color='b', kind='efficiency_units', log=False, 
                              reset=True, **fig_kw):
        """
        Plots an impulse response function.
        
        Arguments:
            
            variables: (list) List of variables whose impulse response functions
                       you wish to plot. Alternatively, you can plot irfs for 
                       all variables by setting variables = 'all'.
                
            method:    (str) Method used to compute the irfs. Must be one of 
                       'linearization', 'forward_shooting', 'reverse_shooting'
                       or 'multiple_shooting.'
                    
            param:     (string) Model parameter.
            
            shock:     (float) Shock to the parameter. Values < 1 correspond to 
                       a reduction in the current value of the parameter; 
                       values > 1 correspond to an increase in the current value
                       of the parameter.
                   
            T:         (float) Length of the impulse response. Default is 40.
            
            kind:      (str) Whether you want impulse response functions in 
                       'levels', 'per_capita', or 'efficiency_units'. Default is
                       for irfs to be in 'efficiency_units'.
                   
            log:   (boolean) Whether or not to have logarithmic scales on the
                   vertical axes. Default is False.
                   
            reset: (boolean) Whether or not to reset the original parameters to
                   their pre-shock values. Default is True.
            
        Returns: A list containing...
        
        """
        # first need to generate and irf
        irf = self.get_impulse_response(method, param, shock, T, kind, reset)
        
        # create mapping from variables to column indices
        irf_dict = {'k':irf[:,[0,1]], 'c':irf[:,[0,2]], 'y':irf[:,[0,3]], 
                    'r':irf[:,[0,4]], 'w':irf[:,[0,5]]}
        
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
            self.plot_trajectory(traj, color, axes[i,0], phase_space=False)
                
            # adjust axis limits
            axes[i,0].set_ylim(0.95 * traj[:,1].min(), 1.05 * traj[:,1].max())
            axes[i,0].set_xlim(2000, 2013 + T)
                
            # y axis labels depend on kind of irfs
            if kind == 'per_capita':
                  
                ti = traj[:,0] - self.data.index[0]
                gr = self.params['g']
                    
                if var in ['k', 'y', 'c']:                     
                    axes[i,0].plot(traj[:,0], traj[0,1] * np.exp(gr * ti), 'k--')
                    axes[i,0].set_ylabel(r'$\frac{%s}{L}(t)$' %var.upper(), 
                                         rotation='horizontal', fontsize=15, 
                                         family='serif')
                elif var == 'w': 
                    axes[i,0].plot(traj[:,0], traj[0,1] * np.exp(gr * ti), 'k--')
                    axes[i,0].set_ylabel('$%s(t)$' %var.upper(), 
                                         rotation='horizontal', fontsize=15, 
                                         family='serif')
                elif var == 'r':
                    axes[i,0].set_ylabel('$%s(t)$' %var, 
                                         rotation='horizontal', fontsize=15, 
                                         family='serif')
                                               
            elif kind == 'levels':
                ti = traj[:,0] - self.data.index[0]
                gr = self.params['n'] + self.params['g']
                    
                if var in ['k', 'y', 'c']:   
                    axes[i,0].plot(traj[:,0], traj[0,1] * np.exp(gr * ti), 'k--')
                    axes[i,0].set_ylabel('$%s(t)$' %var.upper(), 
                                         rotation='horizontal', fontsize=15, 
                                         family='serif')
                elif var == 'w':
                    axes[i,0].plot(traj[:,0], traj[0,1] * np.exp(gr * ti), 'k--')
                    axes[i,0].set_ylabel('$%sL(t)$' %var.upper(), 
                                         rotation='horizontal', fontsize=15, 
                                         family='serif')
                elif var == 'r':
                    axes[i,0].set_ylabel('$%s(t)$' %var, 
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
    
    def solve_linearization(self, k0, ti):
        """
        Solves for a first-order Taylor approximation of the stable manifold 
        in a neighborhood of the long-run steady state of the model.
        
        Arguments:
            
            k0: (float) Initial condition for capital per effective worker.
            
            ti: (array) Array of time points at which you wish to compute the 
                linear approximation of k(t) and c(t).
                
        Returns:
            
            linearized_traj: (array) Trajectory representing a linear 
                             approximation of the stable manifold. 
        
        """
        # compute steady state values
        k_star = self.steady_state.values['k_star']
        c_star = self.steady_state.values['c_star']
        
        # initial condition for predetermined var
        ktilde0 = k0 - k_star
         
        # evaluate the jacobian at steady state
        eval_jacobian = self.jac(0, [k_star, c_star], self.params)
        
        # compute eigen*
        eigen_vals, eigen_vecs = linalg.eig(eval_jacobian, left=True, 
                                            right=False)
        
        # get indices of stable/unstable eigenvals  
        stable_idx   = np.where(np.real(eigen_vals) < 0)[0]
        unstable_idx = np.where(np.real(eigen_vals) > 0)[0]
        
        # number of pre-determined variables
        n_y = stable_idx.size
                                                                                
        # decompose matrix of eigenvectors
        P   = eigen_vecs.T
        P10 = P[unstable_idx, :n_y]
        P11 = P[unstable_idx, n_y:]
        
        # decompose eigenvalues
        D0 = eigen_vals[stable_idx]
                        
        # compute the path of the pre-determined variables
        ktilde_traj = ktilde0 * np.exp(D0 * ti[:,np.newaxis])
        
        # compute the path of the free-control variables
        ctilde_traj = -linalg.inv(P11).dot(P10) * ktilde_traj
        
        # undo change of variables
        k_traj = k_star + ktilde_traj
        c_traj = c_star + ctilde_traj
        
        # combine into a model trajectory (cast to real numbers)
        linearized_traj = np.hstack((ti[:,np.newaxis], k_traj, c_traj)).real
        
        return linearized_traj
        
    def solve_forward_shooting(self, k0, h=1e-1, tol=1.5e-4, mesg=False, 
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
            c_h = self.__get_k_locus(0, [k0], self.params)[0]
        else:
            c_l = self.__get_k_locus(0, [k0], self.params)[0]
            c_h = ((1 - self.params['delta']) * k0 + 
                   self.output(0, k0, self.params))

        # default initial guess for c 
        c0 = (c_h + c_l) / 2
        
        # create an instance of the scipy.integrate.ode class      
        ode  = integrate.ode(self.f, self.jac)
        
        # select the integrator
        ode.set_integrator(integrator, **kwargs)
        
        # pass the model parameters as additional args
        ode.set_f_params(*self.args)
        ode.set_jac_params(*self.args)
        
        # set the initial condition
        ode.set_initial_value([k0, c0], 0)
        
        ########## Compute the optimal c0 using bisection method ############
        while ode.successful():
            # integrate the system one step
            ode.integrate(ode.t + h)
            
            # get the values of the vector field
            k_dot, c_dot = self.f(ode.t, ode.y, *self.args)
            
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
        
    def solve_reverse_shooting(self, k0, h=1e-3, eps=1e-5, integrator='lsoda', 
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

        # evaluate the jacobian at steady state
        evaluated_jacobian = self.jac(0, [k_star, c_star], self.params)
        
        # compute eigenvalues and eigenvectors for evaluated_jacobian
        vals, vecs = linalg.eig(evaluated_jacobian)
        
        # find index of the stable eigenvalue
        stable_idx = np.where(np.real(vals) < 0)[0][0]
        
        # local slope of optimal policy evaluated at steady state
        Ck_prime = vecs[1, stable_idx] / vecs[0, stable_idx]
        
        # RHS of equation 10.7.5 from Judd (1998)
        c_prime = lambda k, c, params: (self.c_dot(0, [k, c], params) / 
                                        self.k_dot(0, [k, c], params))
         
        # create an instance of the scipy.integrate.ode class      
        ode = integrate.ode(c_prime)
        
        # select the integrator
        ode.set_integrator(integrator, **kwargs)
        
        # pass the model parameters as additional args
        ode.set_f_params(*self.args)
        
        if k0 > k_star:
            # set initial conditions
            ode.set_initial_value(c_star + eps * Ck_prime, k_star + eps)
            
            # create a storage container for the trajectory
            solution = np.hstack((ode.t, ode.y)) 
        
            # generate a solution trajectory
            while ode.successful():               
                ode.integrate(ode.t + h, step, relax)
                current_step = np.hstack((ode.t, ode.y))
                solution = np.vstack((solution, current_step))  
            
                # check to see if initial condition has been reached
                if ode.t > k0:
                    break
                else:
                    continue 
        else:
            # set initial condition
            ode.set_initial_value(c_star - eps * Ck_prime, k_star - eps)
            
            # create a storage container for the trajectory
            solution = np.hstack((ode.t, ode.y)) 
        
            # generate a solution trajectory
            while ode.successful():               
                ode.integrate(ode.t - h, step, relax)
                current_step = np.hstack((ode.t, ode.y))
                solution = np.vstack((solution, current_step))  
            
                # check to see if initial condition has been reached
                if ode.t < k0:
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
        f   = lambda t, vec: self.f(t, vec, *self.args)
        jac = lambda t, vec: self.jac(t, vec, *self.args)
                
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
        Computes the stable manifold for the optimal growth model.
        
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
        if method == 'reverse_shooting':
            
            # extract the keyword args for 'reverse' shooting
            h           = kwargs.get('h', 1.0)
            eps         = kwargs.get('eps', 1e-6)
            integrator  = kwargs.get('integrator', 'dopri')
            options     = kwargs.get('options', {})
            
            lower_M_S = self.solve_reverse_shooting(kmin, h, eps, integrator, 
                                                    **options)
            upper_M_S = self.solve_reverse_shooting(kmax, h, eps, integrator, 
                                                    **options)
          
            # reverse the direction of lower_M_S
            lower_M_S = lower_M_S[::-1]
        
        elif method == 'forward_shooting':
            
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
            
        elif method == 'multiple_shooting':
            
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

    def get_consumption_policy(self, kmin=None, kmax=None, method='reverse', 
                               deg=3, **kwargs):
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
                                 
            method:     (str) One of 'linearization', 'forward_shooting', 
                        'reverse_shooting', or 'multiple_shooting', depending.
                        
            **kwargs:   (dict) Dictionary of method specific keyword arguments. 
    
                        For method = 'forward_shooting' the following keyword 
                        arguments are required:
                            
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
            
                            boundary_conditions: (callable): A function which 
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
        if method in ['forward_shooting', 'reverse_shooting', 'multiple_shooting']:    
        
            # compute the stable manifold
            M_S = self.get_stable_manifold(kmin, kmax, method, **kwargs)
        
            # construct a callable B-spline repr of the policy function
            c_k = interpolate.UnivariateSpline(M_S[:,0], M_S[:,1], k=deg, s=0)                   
        
        elif method == 'linearization':
            
            # compute steady state values
            k_star = self.steady_state.values['k_star']
            c_star = self.steady_state.values['c_star']
                 
            # evaluate the jacobian at steady state
            eval_jacobian = self.jac(0, [k_star, c_star], self.params)
        
            # compute eigen*
            eigen_vals, eigen_vecs = linalg.eig(eval_jacobian, left=True, 
                                                right=False)
        
            # get indices of stable/unstable eigenvals  
            stable_idx   = np.where(np.real(eigen_vals) < 0)[0]
            unstable_idx = np.where(np.real(eigen_vals) > 0)[0]
        
            # number of pre-determined variables
            n_y = stable_idx.size
                                                                                
            # decompose matrix of eigenvectors
            P   = eigen_vecs.T
            P10 = P[unstable_idx, :n_y]
            P11 = P[unstable_idx, n_y:]                        
        
            # compute the consumption policy function
            c_k = lambda k: c_star - linalg.inv(P11).dot(P10) * (k - k_star)
     
        else:
            raise ValueError, "Invalid 'method'!" 
                              
        return c_k

    def compare_policies(self, pol1, pol2, grid, metric='L2'):
        """
        Compare two consumption policy functions at common grid points.
        
        Arguments:
            
            pol1:   (callable) Consumption policy function computed using the 
                    get_consumption_policy method.
                  
            pol2:   (callable) Consumption policy function computed using the 
                    get_consumption_policy method.
                  
            grid:   (array-like) Grid of values for capital (per person/
                    effective person) at which to compare the two policy 
                    functions.
                  
            metric: (str) One of 'L2', 'maximal', or None depending on whether 
                    you wish to compute L2 errors, maximal errors, or the simple
                    element-wise difference difference between the two policies.
                  
        Returns:
            
            error: (array) Returns the L2 error, maximal error, or element-wise
                   difference between the two policies.
        
        """
        diff = pol1(grid) - pol2(grid)
        
        if metric == 'L2':
            error = np.sum(diff**2)
        elif metric == 'maximal':
            error = np.abs(diff)
        elif metric == None:
            error = diff
        else:
            ValueError
            
        return error
        
def calibrate_ces(model, iso3_code, bounds=None, sigma0=1.0, alpha0=None, 
                  theta0=2.5, rho=0.04):
    """
    Calibrates an optimal growth model with constant elasticity of substition 
    (CES) production using data from the Penn World Tables (PWT).

    Arguments:
        
        model:     (object) An instance of the Model class.
            
        iso3_code: (str) A valid ISO3 country code.        
                   
        bounds:    (tuple) Start and end years for the subset of the PWT data
                   to use for calibration. Default is None which results in all
                   available data being using in calibration.
                   
        sigma0:    (float) User specified initial guess of the true value for 
                   the elasticity of substitution between capital and effective 
                   labor. Setting sigma0 = 1.0 will calibrate a model with 
                   Cobb-Douglas production. If a value of sigma0 != 1.0 is
                   provided, then sigma (and alpha) will be jointly calibrated 
                   using a non-linear least squares approach.
                   
        alpha0:    (float) User specified initial guess of the true value for
                   alpha. If sigma0 = 1.0, then alpha will be calibrated using 
                   data on capital's share. If sigma0 != 1.0, then alpha0 will
                   be used as an initial guess for the non-linear least squares
                   routine used to jointly calibrate sigma and alpha.
                   
        theta0:    (float) User specified initial guess of the true value for 
                   the inverse elasticity of intertemporal substitution.
                   
    TODO: Validate the non-linear least squares approach being used to find
          optimal alpha and sigma. This function needs a refactoring!
        
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
        
    ##### estimate the labor force growth rate #####
        
    # regress log employed persons on linear time trend
    N     = model.data.index.size
    trend = pd.Series(np.linspace(0, N - 1, N), index=model.data.index)
    res   = pd.ols(y=np.log(model.data.emp.loc[start:end]), 
                   x=trend.loc[start:end])
    n     = res.beta[0]
    L0    = np.exp(res.beta[1])
                           
    ##### estimate the depreciation rate for total capital #####
        
    # use average depreciation rate for total capital
    delta = model.dep_rates.delta_k.loc[start+1:end].mean()  
       
    ##### estimate alpha, sigma, and g #####
    
    # extract some data series that will be needed
    A = model.data.rtfpna
    Y = model.data.rgdpna
    K = model.data.rkna
    L = model.data.emp
    H = model.data.hc
            
    if sigma0 == 1.0:
        
        # Cobb-Douglas production
        sigma = 1.0
        
        # with Cobb-Douglas production alpha is simply capital's share
        alpha = (1 - model.data.labsh.loc[start:end]).mean() 
                
        # adjust measure of TFP to remove human capital contribution
        model.data['atfpna'] = A**(1 / (1 - alpha)) * H
        
        # regress log TFP on linear time trend
        res   = pd.ols(y=np.log(model.data.atfpna.loc[start:end]), x=trend)
        g     = res.beta[0]
        A0    = np.exp(res.beta[1])
        
    else:
        
        def objective(params):
            """Optimize to find alpha and sigma."""
            alpha = params[0]
    	    sigma = params[1]
            gamma = (sigma - 1) / sigma
    
            A = model.data.rtfpna.loc[start:end]
            Y = model.data.rgdpna.loc[start:end]
            K = model.data.rkna.loc[start:end]
            L = model.data.emp.loc[start:end]
            H = model.data.hc.loc[start:end]
    
            # adjust measure of TFP to remove human capital contribution 
            A_tilde = ((K / L) * (((A**gamma * (K / L)**(-gamma * (1 - alpha)) * 
                       H**(gamma * (1 - alpha))) - alpha) / (1 - alpha))**(1 / gamma)) 
    
            y = (Y / (A_tilde * L))
            k = (K / (A_tilde * L))
    
            # remove the discontinuity 
            if abs(gamma) < 1e-6:
                rss = np.log(y) - alpha * np.log(k)
            else:
                rss = np.log(y) - (1 / gamma) * np.log(alpha * k**gamma + (1 - alpha))
            return rss
            
        # solve the non-linear least squares problem         
        res = optimize.leastsq(objective, [alpha0, sigma0], full_output=True)
        
        # extract parameter estimates
        alpha, sigma = res[0]
        
        # print the total sum of squares
        tss = np.sum(res[2]['fvec']**2)
        print 'Total sum of squares:', tss
        
        # print info about estimation success/failure
        if res[4] in [1, 2, 3, 4]:
            print 'Sucessful estimation of sigma and alpha:', res[3]
        else:
            print 'Estimation of sigma and alpha failed!', res[3]      
            print 'Parameters will be calibrated using your initial condition.'
             
        # adjust measure of TFP to remove human capital contribution    
        gamma = (sigma - 1) / sigma
        model.data['atfpna'] = ((K / L) * (((A**gamma * (K / L)**(-gamma * (1 - alpha)) * 
                                H**(gamma * (1 - alpha))) - alpha) / (1 - alpha))**(1 / gamma)) 
        
        # regress log TFP on linear time trend
        res   = pd.ols(y=np.log(model.data.atfpna.loc[start:end]), 
                       x=trend.loc[start:end])
                   
        g     = res.beta[0]
        A0    = np.exp(res.beta[1])
    
    #### estimate the coefficient or relative risk aversion #####
    
    # theta is chosen in order to match average K/Y ratio
    capital_output_ratio = (model.data.rkna.loc[start:end] / 
                            model.data.rgdpna.loc[start:end])
    target1 = capital_output_ratio.mean()
        
    def objective2(theta):
        """Defines model predicted capital-output ratio on BGP."""
        # define a dictionary of temporary params
        tmp_params = {'alpha':alpha, 'delta':delta, 'rho':rho, 'n':n, 'g':g,
                      'theta':theta, 'sigma':sigma}
                      
        # extract tmp bgp values
        tmp_k_star = model.steady_state.functions['k_star'](tmp_params)
        tmp_y_star = model.output(0, tmp_k_star, tmp_params)

        # compare model prediction to target
        out = (tmp_k_star / tmp_y_star) - target1
            
        return out
    
    res   = optimize.root(objective2, theta0, method='hybr')
    
    # extract the parameter results
    theta = res.x[0]
    
    if res.success == True:
        print 'Calibration of theta successful!', res.message
    else:
        print 'Calibration of theta failed:', res.message
        print 'Parameter will be calibrated using your initial condition.'
         
    # create a dictionary of model parameters
    params = {'alpha':alpha, 'sigma':sigma, 'delta':delta, 'theta':theta, 
              'rho':rho, 'n':n, 'L0':L0, 'g':g, 'A0':A0}
              
    # update the model's parameters
    model.update_model_parameters(params)
                    
    # compute new steady state values
    model.steady_state.set_values()
