"""
Filename: _model_individual.py

Authors: Lorenzo Lagos

Solves the infinite horizon dynamic learning model 
with value funtion iteration for a single worker.

"""
from textwrap import dedent
from scipy.stats import norm as norm_distribution
from scipy import spatial
import numpy as np
from div0 import div0
from scipy import interp
import scipy.integrate as integrate

class FireProblem(object):
    
    """

    Global Parameters
    ----------
    beta : scalar(float), optional(default=0.95)
        The discount parameter
    b : scalar(float), optional(default=0.358)
        The workers' mandated benefits
    f : scalar(float), optional(default=0.04)
        The firing fines imposed on firms

    Grid Parameters
    ----------
    T_max : scalar(int), optional(default=360)
        The number of periods in the tenure grid (divide by 10 to get months)
    sp : scalar(float), optional(default=2)
        The multiples of std. dev. on the grid of firm beliefs
    st : scalar(float), optional(default=0.1)
        The step on the grid of firm beliefs

    Occupation Parameters
    ----------
    T_1 : scalar(int), optional(default=0)
        The end of the probationary period one
    T_2 : scalar(int), optional(default=30)
        The end of the probationary period two
    y_0 : scalar(float), optional(default=-1)
        The difference between prior expected match quality and wages
    sig2_0 : scalar(float), optional(default=0.5)
        The prior on match quality variance
    sig2_star : scalar(float), optional(default=2)
        The noise of the match quality signal

    Match Variables
    ----------
    w : scalar(float), optional(default=40)
        The monthly wage of the workers
    H_0 : scipy.stats._distn_infrastructure.rv_frozen
        The prior/objective distribution of match quality
    h_0 : function
        The prior/objective density of match quality
    y_star : np.ndarray
        Realized match quality

    Grid Setup
    ----------
    n : scalar(int)
        The number of workers being simulated (size of w)
    y : np.ndarray
        The grid of firms' belief about match quality
    t : np.ndarray
        The grid of worker tenure
    mu_grid : np.ndarray
        The mean belief about future match quality conditional on current beliefs, ndim = n x y.size x t.size
    sigma_grid : np.ndarray
        The variance in beliefs about future match quality conditional on current beliefs, ndim = n x y.size x t.size
    h_grid : np.ndarray
        The density of belief about future match quality conditional on current beliefs, ndim = n x y.size x t.size


    Tenure Variables
    ----------
    c : np.ndarray
        Firing costs, ndim = n x T_max
    costs : np.ndarray
        Firing costs, ndim = n x y.size x t.size

    """

    def __init__(self, w=40, beta=0.95, y_0=-1, sig2_0=0.5, 
                 sig2_star=2, T_max=360, T_1=0, T_2=30,
                 b=0.358, f=0.04, sp=2, st=0.1):
        "Initializing match relations with probationary contracts"
        

        #Initialize
        self.w = w
        self.beta, self.T_max, self.T_1, self.T_2 = beta, T_max, T_1, T_2
        self.y_0, self.sig2_0, self.sig2_star = y_0, sig2_0, sig2_star
        self.b, self.f, self.sp, self.st = b, f, sp, st

        #Match level variables
        self.H_0 = norm_distribution(y_0+(w/10), sig2_0**(0.5))
        self.h_0 = self.H_0.pdf
        self.y_star = self.H_0.rvs(1)

        #Grid setup
        self.n = 1
        self.y = np.arange((y_0+w/10)-(np.ceil(100*(sig2_0)**(0.5))/100)*sp, 
             (y_0+w/10)+(np.ceil(100*(sig2_0)**(0.5))/100)*sp, st) # Belief grid
        self.t = np.arange(1, T_max+1, 1) # Tenure grid starts at 1
        self.mu_grid = np.array([((sig2_0)/(self.t*sig2_0+sig2_star))*self.y_star+
             (((self.t-1)*sig2_0+sig2_star)/(self.t*sig2_0+sig2_star))*k for k in self.y]).reshape(self.n,self.y.size,T_max)
        self.mu_grid = np.c_[np.repeat(y_0+(w/10),self.y.size).reshape(1,self.y.size,1), self.mu_grid] #include priors
        self.sigma_grid = np.tile(((sig2_0)/(self.t*sig2_0+sig2_star))**(0.5)*sig2_star, (self.n,self.y.size,1))
        self.sigma_grid = np.c_[np.repeat(sig2_0,self.y.size).reshape(1,self.y.size,1), self.sigma_grid] #include priors
        self.h_grid = np.array([norm_distribution(i,j**(0.5)).pdf for i, j in zip(self.mu_grid.flatten(),self.sigma_grid.flatten())])

        #Tenure variables
        self.c = np.hstack((np.array(0.5*(T_1-self.t[self.t<=T_1])*(w/10)), 
             np.array(0.5*(T_2-self.t[(self.t>T_1)&(self.t<=T_2)])*(w/10)), 
             np.array(w*(1+b+f*((self.t[self.t>T_2]+10)/10)))))
        self.costs = np.array([np.tile(self.c,(self.y.size,1)) for k in range(self.n)])
        self.costs = np.c_[np.repeat(0,self.y.size).reshape(1,self.y.size,1), self.costs] #include period zero

    def __repr__(self):
        "Meant for the user of an application to see"

        m = "FireProblem(w={wa}, beta={b}, y_0={y0}, sig2_0={s20}, sig2_star={s2s}, "
        m += "T_max={tm}, T_1={t1}, T_2={t2}, b={ben}, f={fin})"
        return m.format(wa=self.w, b=self.beta, y0=self.y_0, s20=self.sig2_0,
                        s2s=self.sig2_star, tm=self.T_max, t1=self.T_1,
                        t2=self.T_2, ben=self.b, fin=self.f)

    def __str__(self):
        "Meant for the programmer to see, as in debugging and development"

        m = """\
        Dynamic learning problem
          - w (monthly wage)                   : {wa}
          - beta (discount parameter)          : {b}
          - y_0 (prior on quality minus wage)  : {y0}
          - H_0 (prior quality)                : Norm(y0+w/10,{s20})
          - sig2_star (noise of signal)        : {s2s}
          - T_max (time horizon in grid)       : {tm}
          - T_1 (length of prob period 1)      : {t1}
          - T_2 (length of prob period 2)      : {t2}
          - b (worker benfits)                 : {ben}
          - f (firing fines)                   : {fin}
          - n (number of workers)              : {num}
        """
        hm, hs = self.H_0.args
        return dedent(m.format(wa=self.w, b=self.beta, y0=self.y_0, s20=self.sig2_0,
                               s2s=self.sig2_star, tm=self.T_max, t1=self.T_1,
                               t2=self.T_2, ben=self.b, fin=self.f, num=self.n))

    def bellman_operator(self, v):
        """
        The Bellman operator for the dynamic learning model.

        Parameters
        ----------
        v : array_like(float)
            A 3D NumPy array representing the value function
            Interpretation: :math:`v[k, i, j] = v(\n_k, \y_i, \t_j)`

        Returns
        -------
        new_v : array_like(float)
            The updated value function Tv as an array of shape v.shape

            """

        new_v = np.empty(v.shape)

        v_cols = np.empty(self.n*self.y.size*(self.t.size+1), dtype=object)
        it=0
        for k in range(self.n):
            v_cols0 = np.array([(lambda z: (lambda x: interp(x, self.y, z)))(new_v[k,:,0])])
            for i in range(self.y.size):
                for j in range(self.t.size+1):
                    if j != self.T_max: #from 0 to (T_max-1)
                        v_cols[it] = np.array([(lambda z: (lambda x: interp(x, self.y, z)))(v[k,:,(j+1)])])
                    else:          #T_max duplicate
                        v_cols[it] = np.array([(lambda z: (lambda x: interp(x, self.y, z)))(v[k,:,j])])
                    it = it+1

        #keep worker
        E0 = np.array([integrate.fixed_quad(lambda x: v_cols[l][0](x)*self.h_grid[l](x), i-3*j**0.5, i+3*j**0.5)[0] 
             for l, (i, j) in enumerate(zip(self.mu_grid.flatten(),self.sigma_grid.flatten()))]).reshape(self.n,self.y.size,self.t.size+1)
        v0 = self.beta*E0

        #Fire worker
        E1 = np.array([integrate.fixed_quad(lambda x: v_cols0[0](x)*self.h_0(x), 
             self.y_0+self.w/10-3*self.sig2_0**0.5, self.y_0+self.w/10+3*self.sig2_0**0.5)[0] 
             for l in range(self.n) for (i, j) in enumerate(zip(self.mu_grid.flatten(),self.sigma_grid.flatten()))]).reshape(self.n,self.y.size,self.t.size+1)
        #E1 = np.array([integrate.fixed_quad(lambda x: v_cols[(self.n*self.y.size*(self.t.size+1))*l][0](x)*self.h_0(x), 
        #     self.y_0+self.w/10-3*self.sig2_0**0.5, self.y_0+self.w/10+3*self.sig2_0**0.5)[0] 
        #     for l in range(self.n) for (i, j) in enumerate(zip(self.mu_grid.flatten(),self.sigma_grid.flatten()))]).reshape(self.n,self.y.size,self.t.size+1)
        v1 = -self.costs + self.beta*E1

        new_v = np.maximum(v0, v1)

        return new_v

    def get_greedy(self, v):
        """
        Compute optimal actions taking v as the value function.

        Parameters
        ----------
        v : array_like(float)
            A 3D NumPy array representing the value function
            Interpretation: :math:`v[k, i, j] = v(n_k, y_i, t_j)`

        Returns
        -------
        policy : array_like(float)
            A 3D NumPy array, where policy[k, i, j] is the optimal action
            at :math:`(n_k, y_i, t_j)`.

            The optimal action is represented as an integer in the set
            [0, 1] where 0 = 'keep' and 1 = 'fire'

        """

        policy = np.empty(v.shape, dtype=int)

        v_cols = np.empty(self.n*self.y.size*(self.t.size+1), dtype=object)
        it=0
        for k in range(self.n):
            v_cols0 = np.array([(lambda z: (lambda x: interp(x, self.y, z)))(new_v[k,:,0])])
            for i in range(self.y.size):
                for j in range(self.t.size+1):
                    if j != self.T_max: #from 0 to (T_max-1)
                        v_cols[it] = np.array([(lambda z: (lambda x: interp(x, self.y, z)))(v[k,:,(j+1)])])
                    else:          #T_max duplicate
                        v_cols[it] = np.array([(lambda z: (lambda x: interp(x, self.y, z)))(v[k,:,j])])
                    it = it+1

        #keep worker
        E0 = np.array([integrate.fixed_quad(lambda x: v_cols[l][0](x)*self.h_grid[l](x), i-3*j**0.5, i+3*j**0.5)[0] 
             for l, (i, j) in enumerate(zip(self.mu_grid.flatten(),self.sigma_grid.flatten()))]).reshape(self.n,self.y.size,self.t.size+1)
        v0 = self.beta*E0

        #Fire worker
        E1 = np.array([integrate.fixed_quad(lambda x: v_cols0[0](x)*self.h_0(x), 
             self.y_0+self.w/10-3*self.sig2_0**0.5, self.y_0+self.w/10+3*self.sig2_0**0.5)[0] 
             for l in range(self.n) for (i, j) in enumerate(zip(self.mu_grid.flatten(),self.sigma_grid.flatten()))]).reshape(self.n,self.y.size,self.t.size+1)
        #E1 = np.array([integrate.fixed_quad(lambda x: v_cols[(self.n*self.y.size*(self.t.size+1))*l][0](x)*self.h_0(x), 
        #     self.y_0+self.w/10-3*self.sig2_0**0.5, self.y_0+self.w/10+3*self.sig2_0**0.5)[0] 
        #     for l in range(self.n) for (i, j) in enumerate(zip(self.mu_grid.flatten(),self.sigma_grid.flatten()))]).reshape(self.n,self.y.size,self.t.size+1)
        v1 = -self.costs + self.beta*E1

        policy = (v1>v0).astype(int)

        return policy
