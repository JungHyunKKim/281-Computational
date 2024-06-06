## SeqSpaceDurable.py 

# Import packages
import numpy as np
import scipy as sp
import time 

# Solve steady state of the model using the Newton-Raphson method
class solveSeqSpaceDurable:
    
    # Define the parameters
    def __init__(self, 
        sigma    = 1, 
        gamma    = 1, 
        beta     = 1/(1+0.01), 
        delta    = 0.017, 
        eta      = 1, 
        DC_ratio = 9.058,
        horizon  = 100):

        self.sigma    = sigma
        self.gamma    = gamma
        self.beta     = beta
        self.delta    = delta
        self.eta      = eta
        self.DC_ratio = DC_ratio
        self.horizon  = horizon

    def solve_steady_state(self): 
        self.rss = 1/self.beta - 1
        self.Qss = 1 
        self.pxss = 1 
        self.Dss = 1/(1/self.DC_ratio + self.delta)
        self.Css = 1 - self.delta*self.Dss
        self.Xss = self.delta*self.Dss
        self.alpha = self.Css**(-self.sigma) * self.Dss**(self.gamma) * (1 - self.beta*(1-self.delta))
        self.Yxss = self.Xss
        self.Nxss = self.Yxss 
        self.Ycss = self.Css
        self.Ncss = self.Ycss
        self.Zss  = 1
        self.wss  = 1

    # Solve the model using the sequential space method
    def solve_seqspace(self, shocks): 

        start = time.time()
        
        # Order of variables: Nc, Yx, Yc, w, X, C, r, D, Q 

        # define sparse identity, above-diagonal sparse matrix, below-diagonal sparse matrix, and zero matrix
        I   = sp.sparse.eye(self.horizon)
        Ip1 = sp.sparse.diags([np.ones(self.horizon-1)], [1], (self.horizon, self.horizon))
        Im1 = sp.sparse.diags([np.ones(self.horizon-1)], [-1], (self.horizon, self.horizon))
        Z   = sp.sparse.csr_matrix((self.horizon, self.horizon)) # sparse empty matrix of size self.horizon x self.horizon
        I0 = sp.sparse.diags([np.hstack(([1], np.zeros(self.horizon-1)))], [0], (self.horizon, self.horizon))

        # Durable Euler equation block
        PhiDEE_Nc = Z
        PhiDEE_Yx = Z
        PhiDEE_Yc = Z
        PhiDEE_w  = Z
        PhiDEE_X  = Z
        PhiDEE_C  = self.sigma*(1-self.beta*(1-self.delta))*I
        PhiDEE_r  = -self.beta*(1-self.delta)*I
        PhiDEE_D  = -self.gamma*(1-self.beta*(1-self.delta))*I
        PhiDEE_Q  = self.beta*Ip1 - I

        # combine matrix blocks in a single sparse matrix with the following structure:
        dHdY = sp.sparse.bmat([[PhiDEE_Nc, PhiDEE_Yx, PhiDEE_Yc, PhiDEE_w, PhiDEE_X, PhiDEE_C, PhiDEE_r, PhiDEE_D, PhiDEE_Q]])

        assert dHdY.shape == (1*self.horizon, 9*self.horizon)

        # Construct dY/dU and dY/dZ using the blocks. 

        ### Labor Market Clearing Block ###
        # Labor market clearing: non-durable sector employment
        Phi_NcNx = -self.delta*self.Dss / (1-self.delta*self.Dss) * I 
        Phi_NcZ  = Z
        Phi_NcD0 = Z

        # combine labor market clearing block matrices
        dYLMCdU = sp.sparse.bmat([[Phi_NcNx]])
        dYLMCdZ = sp.sparse.bmat([[Phi_NcZ, Phi_NcD0]])

        assert dYLMCdU.shape == (1*self.horizon, 1*self.horizon)
        assert dYLMCdZ.shape == (1*self.horizon, 2*self.horizon)

        ### Firm block matrices ###

        ## firm block matrices: durable output
    
        # total effect
        Phi_YxNx = I
        Phi_YxZ  = I  
        Phi_YxD0 = Z 

        ## firm block matrices: non-durable output
        
        # intermediate effect
        PhiIE_YcNc = I 

        # total effect
        Phi_YcNx = PhiIE_YcNc * Phi_NcNx
        Phi_YcZ  = I + PhiIE_YcNc * Phi_NcZ
        Phi_YcD0 = PhiIE_YcNc * Phi_NcD0

        ## firm block matrices: wage
        
        # total effect
        Phi_wNx = Z
        Phi_wZ = I 
        Phi_wD0 = Z

        # combine firm block matrices
        dYFdU = sp.sparse.bmat([[Phi_YxNx],
                                [Phi_YcNx],
                                [Phi_wNx]])
        
        dYFdZ = sp.sparse.bmat([[Phi_YxZ, Phi_YxD0],
                                [Phi_YcZ, Phi_YcD0],
                                [Phi_wZ, Phi_wD0]])

        assert dYFdU.shape == (3*self.horizon, 1*self.horizon)
        assert dYFdZ.shape == (3*self.horizon, 2*self.horizon)

        ### Goods market clearing ###

        ## goods market clearing: durable expenditure
        # intermediate effect
        PhiIE_XYx = I

        # total effect
        Phi_XNx = PhiIE_XYx * Phi_YxNx
        Phi_XZ  = PhiIE_XYx * Phi_YxZ
        Phi_XD0 = PhiIE_XYx * Phi_YxD0

        ## goods market clearing: non-durable goods

        # intermediate effect
        PhiIE_CYc = I

        # total effect
        Phi_CNx = PhiIE_CYc * Phi_YcNx
        Phi_CZ  = PhiIE_CYc * Phi_YcZ
        Phi_CD0 = PhiIE_CYc * Phi_YcD0

        # combine goods market clearing matrices
        dYGMdU = sp.sparse.bmat([[Phi_XNx],
                                 [Phi_CNx]])
        
        dYGMdZ = sp.sparse.bmat([[Phi_XZ, Phi_XD0],
                                 [Phi_CZ, Phi_CD0]])

        assert dYGMdU.shape == (2*self.horizon, 1*self.horizon)
        assert dYGMdZ.shape == (2*self.horizon, 2*self.horizon)

        ### Household Block ###

        ## household block: interest rate
        # intermediate effect
        PhiIE_rC = self.sigma * (Ip1 - I) 

        # total effect
        Phi_rNx = PhiIE_rC * Phi_CNx
        Phi_rZ  = PhiIE_rC * Phi_CZ
        Phi_rD0 = PhiIE_rC * Phi_CD0

        ## household block: durable stock
        # intermediate effect
        temp = (I - (1-self.delta)*Im1)
        PhiIE_DX = sp.sparse.linalg.spsolve(temp, self.delta*I)
        PhiIE_DD0 = sp.sparse.linalg.spsolve(temp, (1-self.delta)*I0)
        
        # total effect
        Phi_DNx = PhiIE_DX * Phi_XNx
        Phi_DZ  = PhiIE_DX * Phi_XZ
        Phi_DD0 = PhiIE_DX * Phi_XD0 + PhiIE_DD0 

        ## household block: real shawdow price of durable stock
        # intermediate effect
        PhiIE_QX = self.eta * I 
        PhiIE_QD = -self.eta * Im1 
        PhiIE_QD0 = -self.eta * I0

        # total effect
        Phi_QNx = PhiIE_QX * Phi_XNx + PhiIE_QD * Phi_DNx
        Phi_QZ  = PhiIE_QX * Phi_XZ + PhiIE_QD * Phi_DZ
        Phi_QD0 = PhiIE_QX * Phi_XD0 + PhiIE_QD * Phi_DD0 + PhiIE_QD0

        # combine household block matrices
        dYHHdU = sp.sparse.bmat([[Phi_rNx],
                                 [Phi_DNx],
                                 [Phi_QNx]])

        dYHHdZ = sp.sparse.bmat([[Phi_rZ, Phi_rD0],
                                 [Phi_DZ, Phi_DD0],
                                 [Phi_QZ, Phi_QD0]])

        assert dYHHdU.shape == (3*self.horizon, 1*self.horizon)
        assert dYHHdZ.shape == (3*self.horizon, 2*self.horizon)

        ### Combine all blocks ###

        # stack dYLMCdU, dYFdU, dYGMdU, dYHHdU
        dYdU = sp.sparse.bmat([[dYLMCdU],
                               [dYFdU],
                               [dYGMdU],
                               [dYHHdU]])
        
        # stack dYLMCdZ, dYFdZ, dYGMdZ, dYHHdZ
        dYdZ = sp.sparse.bmat([[dYLMCdZ],
                               [dYFdZ],
                               [dYGMdZ],
                               [dYHHdZ]])
        
        assert dYdU.shape == (9*self.horizon, 1*self.horizon)
        assert dYdZ.shape == (9*self.horizon, 2*self.horizon)

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
        assert dYdZ.shape == (9*self.horizon, 2*self.horizon)
        assert dXdZ.shape == (10*self.horizon, 2*self.horizon)

        # compute impulse response functions
        X = dXdZ @ shocks

        print(f'Elapsed time is {time.time()-start:.2f} seconds.')    

        return X

