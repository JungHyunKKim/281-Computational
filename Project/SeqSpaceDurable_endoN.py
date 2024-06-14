## SeqSpaceDurable_endoN.py 

# Import packages
import numpy as np
import scipy as sp
import time 

class solveSeqSpaceDurable:
    
    # Define the parameters
    def __init__(self, 
        sigma    = 1, 
        gamma    = 1, 
        beta     = 1/(1+0.01), 
        delta    = 0.017, 
        eta      = 1, 
        DC_ratio = 9.058,
        chi      = 1, 
        psi      = 1, 
        horizon  = 100):

        self.sigma    = sigma
        self.gamma    = gamma
        self.beta     = beta
        self.delta    = delta
        self.eta      = eta
        self.DC_ratio = DC_ratio
        self.chi      = chi
        self.psi      = psi
        self.horizon  = horizon

    # Define the excess non-durable labor demand function
    def excess_Nc(self, Ncss, Dss):

        # Excess non-durable labor demand in SS : Nc = N - Nx
        excess_Nc = 1/(self.chi**(1/self.psi) * Ncss**(self.sigma/self.psi)) - self.delta * Dss - Ncss

        return excess_Nc

    # Derivative
    def excess_Nc_derivative(self, Ncss, Dss, epsilon=1e-7): 
        return (self.excess_Nc(Ncss, Dss) - self.excess_Nc(Ncss - epsilon, Dss)) / (epsilon)

    # Newton-Raphson method
    def NR(self, f, df, x0, tol=1e-8, max_iter_NR=1000):
        x = x0
        for i in range(max_iter_NR):
            x_new = x - f(x) / df(x)
            if abs(x_new - x) < tol:
                break
            x = x_new
        return x_new

    def solve_steady_state(self): 
        self.rss = 1/self.beta - 1
        self.Qss = 1 
        self.pxss = 1 

        if self.chi == 0:
            self.Dss = 1/(1/self.DC_ratio + self.delta)
        else:
            self.Dss = 1/( self.chi**(1/self.psi) * (1/self.DC_ratio)**(self.sigma/self.psi) * (1/self.DC_ratio + self.delta))**(self.psi/(self.sigma + self.psi))
        
        self.Xss = self.delta*self.Dss
        self.Yxss = self.Xss
        self.Nxss = self.Yxss 
        self.Zss  = 1
        self.wss  = 1

        # define excess non-durable labor demand function evaluated at self.Dss
        if self.chi == 0:
            self.Ncss = 1 - self.delta*self.Dss
        else:
            self.excess_Nc_temp = lambda Ncss: self.excess_Nc(Ncss, self.Dss)
            self.excess_Nc_derivative_temp = lambda Ncss: self.excess_Nc_derivative(Ncss, self.Dss)
            self.Ncss = self.NR(self.excess_Nc_temp, self.excess_Nc_derivative_temp, 0.2)
        
        self.Css = self.Ncss
        self.Ycss = self.Ncss

        if self.chi ==0:
            self.Nss = 1 
        else: 
            self.Nss = 1/(self.chi**(1/self.psi) * self.Ncss**(self.sigma/self.psi))

        self.alpha = self.Css**(-self.sigma) * self.Dss**(self.gamma) * (1 - self.beta*(1-self.delta))
        self.omega_c = self.Ncss / self.Nss
        self.omega_x = 1 - self.omega_c

    # Solve the model using the sequential space method
    def solve_seqspace(self, shocks): 

        start = time.time()
        
        # Order of variables: U = Nc, Y = Yc, w, C, r, N, Nx, Yx, X, D, Q

        # define sparse identity, above-diagonal sparse matrix, below-diagonal sparse matrix, and zero matrix
        I   = sp.sparse.eye(self.horizon)
        Ip1 = sp.sparse.diags([np.ones(self.horizon-1)], [1], (self.horizon, self.horizon))
        Im1 = sp.sparse.diags([np.ones(self.horizon-1)], [-1], (self.horizon, self.horizon))
        Z   = sp.sparse.csr_matrix((self.horizon, self.horizon)) # sparse empty matrix of size self.horizon x self.horizon
        I0 = sp.sparse.diags([np.hstack(([1], np.zeros(self.horizon-1)))], [0], (self.horizon, self.horizon))

        # Durable Euler equation block
        PhiDEE_Yc = Z
        PhiDEE_w  = Z
        PhiDEE_C  = self.sigma*(1-self.beta*(1-self.delta))*I
        PhiDEE_r  = -self.beta*(1-self.delta)*I
        PhiDEE_N = Z
        PhiDEE_Nx = Z 
        PhiDEE_Yx = Z
        PhiDEE_X  = Z
        PhiDEE_D  = -self.gamma*(1-self.beta*(1-self.delta))*I
        PhiDEE_Q  = self.beta*Ip1 - I

        # combine matrix blocks in a single sparse matrix with the following structure:
        dHdY = sp.sparse.bmat([[PhiDEE_Yc, PhiDEE_w, PhiDEE_C, PhiDEE_r, PhiDEE_N, PhiDEE_Nx, PhiDEE_Yx, PhiDEE_X, PhiDEE_D, PhiDEE_Q]])

        assert dHdY.shape == (1*self.horizon, 10*self.horizon)

        # Construct dY/dU and dY/dZ using the blocks. 

        ### Non-durable sector block ###

        ## non-durable output
        Phi_YcNc = I
        Phi_YcZ  = I
        Phi_YcD0 = Z

        # wage
        Phi_wNc = Z
        Phi_wZ  = I
        Phi_wD0 = Z

        ## non-durable consumption
        # intermediate effect
        PhiIE_CYc = I 

        # total effect
        Phi_CNc = PhiIE_CYc @ Phi_YcNc
        Phi_CZ  = PhiIE_CYc @ Phi_YcZ
        Phi_CD0 = PhiIE_CYc @ Phi_YcD0

        ## real interest rate
        # intermediate effect
        PhiIE_rC = self.sigma * (Ip1 - I)

        # total effect
        Phi_rNc = PhiIE_rC @ Phi_CNc
        Phi_rZ  = PhiIE_rC @ Phi_CZ
        Phi_rD0 = PhiIE_rC @ Phi_CD0

        ## aggregate employment
        # intermediate effect
        PhiIE_Nw = 1/self.psi * I 
        PhiIE_NC = - self.sigma/self.psi * I

        # total effect
        Phi_NNc = PhiIE_Nw @ Phi_wNc + PhiIE_NC @ Phi_CNc
        Phi_NZ  = PhiIE_Nw @ Phi_wZ + PhiIE_NC @ Phi_CZ
        Phi_ND0 = PhiIE_Nw @ Phi_wD0 + PhiIE_NC @ Phi_CD0

        if self.chi == 0: 
            Phi_NNc = Z
            Phi_NZ  = Z
            Phi_ND0 = Z

        # combine non-durable sector matrices
        dYNDdU = sp.sparse.bmat([[Phi_YcNc],
                                 [Phi_wNc],
                                 [Phi_CNc],
                                 [Phi_rNc],
                                 [Phi_NNc]])

        dYNDdZ = sp.sparse.bmat([[Phi_YcZ, Phi_YcD0],
                                 [Phi_wZ, Phi_wD0],
                                 [Phi_CZ, Phi_CD0],
                                 [Phi_rZ, Phi_rD0],
                                 [Phi_NZ, Phi_ND0]])
        

        ### Labor Market Clearing Block ###

        ## Durable sector employment
        # intermediate effect
        PhiIE_NxN = 1/self.omega_x * I
        PhiIE_NxNc = - self.omega_c/self.omega_x * I

        # total effect
        Phi_NxNc = PhiIE_NxN @ Phi_NNc + PhiIE_NxNc
        Phi_NxZ  = PhiIE_NxN @ Phi_NZ
        Phi_NxD0 = PhiIE_NxN @ Phi_ND0

        # combine labor market clearing block matrices
        dYLMCdU = sp.sparse.bmat([[Phi_NxNc]])
        dYLMCdZ = sp.sparse.bmat([[Phi_NxZ, Phi_NxD0]])

        assert dYLMCdU.shape == (1*self.horizon, 1*self.horizon)
        assert dYLMCdZ.shape == (1*self.horizon, 2*self.horizon)

        ### Durable sector block matrices ###

        ## durable output
        # intermediate effect
        PhiIE_YxNx = I
        PhiIE_YxZ  = I

        # total effect
        Phi_YxNc = PhiIE_YxNx @ Phi_NxNc
        Phi_YxZ  = PhiIE_YxNx @ Phi_NxZ + PhiIE_YxZ
        Phi_YxD0 = PhiIE_YxNx @ Phi_NxD0

        ## durable expenditure
        # intermediate effect
        PhiIE_XYx = I

        # total effect
        Phi_XNc = PhiIE_XYx @ Phi_YxNc
        Phi_XZ  = PhiIE_XYx @ Phi_YxZ
        Phi_XD0 = PhiIE_XYx @ Phi_YxD0

        ## durable stock
        # intermediate effect
        temp = (I - (1-self.delta)*Im1)
        PhiIE_DX = sp.sparse.linalg.spsolve(temp, self.delta*I)
        PhiIE_DD0 = sp.sparse.linalg.spsolve(temp, (1-self.delta)*I0)

        # total effect
        Phi_DNc = PhiIE_DX @ Phi_XNc
        Phi_DZ  = PhiIE_DX @ Phi_XZ
        Phi_DD0 = PhiIE_DX @ Phi_XD0 + PhiIE_DD0

        ## real shadow price of durable stock
        # intermediate effect
        PhiIE_QX = self.eta * I
        PhiIE_QD = -self.eta * Im1
        PhiIE_QD0 = -self.eta * I0

        # total effect
        Phi_QNc = PhiIE_QX @ Phi_XNc + PhiIE_QD @ Phi_DNc
        Phi_QZ  = PhiIE_QX @ Phi_XZ + PhiIE_QD @ Phi_DZ
        Phi_QD0 = PhiIE_QX @ Phi_XD0 + PhiIE_QD @ Phi_DD0 + PhiIE_QD0

        # combine durable sector matrices
        dYDdU = sp.sparse.bmat([[Phi_YxNc],
                                [Phi_XNc],
                                [Phi_DNc],
                                [Phi_QNc]])
        
        dYDdZ = sp.sparse.bmat([[Phi_YxZ, Phi_YxD0],
                                [Phi_XZ, Phi_XD0],
                                [Phi_DZ, Phi_DD0],
                                [Phi_QZ, Phi_QD0]])

        assert dYDdU.shape == (4*self.horizon, 1*self.horizon)
        assert dYDdZ.shape == (4*self.horizon, 2*self.horizon)

        ### Combine all blocks ###

        # stack dYNDdU, dYLMCdU, dYDdU to get dYdU
        dYdU = sp.sparse.bmat([[dYNDdU],
                               [dYLMCdU],
                               [dYDdU]])
        
        # stack dYNDdZ, dYLMCdZ, dYDdZ to get dYdZ
        dYdZ = sp.sparse.bmat([[dYNDdZ],
                               [dYLMCdZ],
                               [dYDdZ]])
        
        assert dYdU.shape == (10*self.horizon, 1*self.horizon)
        assert dYdZ.shape == (10*self.horizon, 2*self.horizon)

        # compute dHdU using the chain rule dHdU = dHdY @ dYdU (@ is the python matrix multiplication operator)
        dHdU = dHdY @ dYdU 

        # compute dHdZ using the chain rule dHdZ = dHdY @ dYdZ (@ is the python matrix multiplication operator)
        dHdZ = dHdY @ dYdZ

        assert sp.sparse.issparse(dHdZ) == True
        assert sp.sparse.issparse(dHdU) == True

        assert dHdU.shape == (1*self.horizon, 1*self.horizon)
        assert dHdZ.shape == (1*self.horizon, 2*self.horizon)

        # compute the Jacobian of the model
        dUdZ = - sp.sparse.linalg.spsolve(dHdU, dHdZ)
        dYdZ = dYdU @ dUdZ + dYdZ

        dXdZ = sp.sparse.bmat([[dUdZ],
                               [dYdZ]])

        assert dUdZ.shape == (1*self.horizon, 2*self.horizon)
        assert dYdZ.shape == (10*self.horizon, 2*self.horizon)
        assert dXdZ.shape == (11*self.horizon, 2*self.horizon)

        # compute impulse response functions
        X = dXdZ @ shocks

        print(f'Elapsed time is {time.time()-start:.2f} seconds.')    

        return X

