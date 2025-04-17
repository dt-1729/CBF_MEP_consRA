import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.linalg import block_diag
from scipy import optimize
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from scipy.spatial.distance import cdist
import testcases
from utils import *
import time
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# create FLP class and its associated functions
class FLP():
    # declare the class variables
    N : int # number of resources
    M : int # number of facilities
    d : int # dimension of the problem
    resLoc : np.ndarray # resource locations
    rho : np.ndarray # rosource weights
    C : np.ndarray # capacity constraints
    B : np.ndarray # min utilization constraints
    P_eps : float # min probability value (> 0)
    beta_tol : float # min beta value (> 0)
    alloc_C : np.ndarray # cost of association
    
    # initialize the class variables
    def __init__(self, N, M, d, resLoc, rho, C, B, alloc_C, P_eps, beta_tol):
        self.N = N
        self.M = M
        self.d = d
        self.resLoc = resLoc
        self.rho = rho
        self.C = C
        self.B = B
        self.P_eps = P_eps
        self.beta_tol = beta_tol
        self.alloc_C = alloc_C
        # print('================= Class FLP() initialized =================')


    ## ================= Clustering DA ================= ##
    
    def get_D(self, Y):
        X = self.resLoc
        D = cdist(X, Y, metric='sqeuclidean')
        return D

    def get_cap_cons(self, P, theta):
        cons = np.sum(self.rho * P * self.alloc_C, axis=0) - self.C
        th_exp = theta * np.exp(theta * cons)
        tiled_th_exp = np.tile(th_exp, (self.N,1))
        return cons, tiled_th_exp

    def get_PY_gibbs(self, D, beta):
        D_min = D.min(axis=1, keepdims=True)
        exp_beta_D = np.exp(-beta * (D-D_min))
        Pi = exp_beta_D/np.sum(exp_beta_D, axis=1, keepdims=True)
        Pi_ex = np.expand_dims(Pi.T, axis=2)
        rho_Pi_X = Pi_ex * (self.resLoc * self.rho)
        rho_Pi = Pi * self.rho
        Yi = np.sum(rho_Pi_X, axis=1)/np.sum(rho_Pi, axis=0).reshape(-1,1)
        return Pi, Yi

    def get_PY_gibbs_C(self, D, P, beta, beta1, theta):
        # D_min = D.min(axis=1, keepdims=True)
        cons, penalty = self.get_cap_cons(P, theta)
        # penalty_min = penalty.min(axis=1, keepdims=True)
        beta_D = beta * D + beta1 * penalty
        beta_D_min = beta_D.min(axis=1,keepdims=True)
        exp_beta_D = np.exp(-(beta_D-beta_D_min))
        sum_exp_beta_D = np.sum(exp_beta_D, axis=1, keepdims=True)
        # print(sum_exp_beta_D)
        Pi = exp_beta_D/sum_exp_beta_D
        Pi_ex = np.expand_dims(Pi.T, axis=2)
        rho_Pi_X = Pi_ex * (self.resLoc * self.rho)
        rho_Pi = Pi * self.rho
        # print(Pi[:,0])
        # print(np.sum(rho_Pi, axis=0).reshape(-1,1))
        Yi = np.sum(rho_Pi_X, axis=1)/np.sum(rho_Pi, axis=0).reshape(-1,1)
        return Pi, Yi

    def get_F(self, Y, beta):
        D = self.get_D(Y)
        P = self.get_PY_gibbs(D, beta)[0]
        F = np.sum(self.rho*np.sum(P*(D + 1/beta * myLog(P)),axis=1,keepdims=True)) + 1/beta * np.log(self.M)
        return F

    def optimize_DA(self, P0, Y0, beta, n_iters, D_tol, allowPrint=False):

        Yi = Y0
        Pi = P0
        for i in range(n_iters):
            D_next = self.get_D(Yi)
            D_next.shape
            Pi, Yi = self.get_PY_gibbs(D_next, beta)
            if i > 0:
                diff_D = D_next - D_prev
                norm_diff_D = np.max(abs(diff_D))
                if allowPrint:
                    print(f'i{i}\tnorm_diff_D:{norm_diff_D:.3e}')
                if norm_diff_D < D_tol:
                    print(f'tolerance achieved: \t norm_diff_D={norm_diff_D:.3e}')
                    break
            D_prev = D_next

        return Pi, Yi
    

    def optimize_DA_C(self, P0, Y0, beta, b1_min, b1_max, b1_grow, theta, n_iters, D_tol, allowPrint1=False, allowPrint2=False):

        Yi = Y0
        Pi = P0
        b1 = b1_min
        while b1 <= b1_max:
            # solve implicit equations
            for i in range(n_iters):
                D_next = self.get_D(Yi)
                Pi, Yi = self.get_PY_gibbs_C(D_next, Pi, beta, b1, theta)
                if i > 0:
                    diff_D = D_next - D_prev
                    norm_diff_D = np.max(abs(diff_D))
                    if allowPrint2:
                        print(f'i{i}\tnorm_diff_D:{norm_diff_D:.3e}')
                    if norm_diff_D < D_tol and allowPrint2:
                        print(f'tolerance achieved: \t norm_diff_D={norm_diff_D:.3e}')
                        break
                D_prev = D_next
            cons_b1 = self.get_cap_cons(Pi, theta)[0]
            b1 = b1 * b1_grow

        cons = self.get_cap_cons(Pi, theta)[0]
        print(f'cons_violation: {np.round(cons,3)}')

        return Pi, Yi
    

    def anneal_DA(self, Y0, b_min, b_max, b_grow, n_iters, D_tol):

        b = b_min
        Yb = Y0
        Db = self.get_D(Yb)
        Pb = self.get_PY_gibbs(Db, b)[0]

        F_array = []
        b_array = []
        P_array = np.array([Pb])
        Y_array = np.array([Yb])
        t_compute_array = []
        Cap_array = np.array([self.get_cap_cons(Pb, theta=1)[0] + self.C])

        while b <= b_max:
            t0 = time.time()
            Pb, Yb = self.optimize_DA(Pb, Yb, b, n_iters, D_tol)
            t1 = time.time()
            Yb += np.random.multivariate_normal(np.zeros(2*self.M), 0.001*np.eye(2*self.M)).reshape(-1,2)
            Fb = self.get_F(Yb, b)
            cap = self.get_cap_cons(Pb, theta=1)[0] + self.C

            F_array.append(Fb)
            b_array.append(b)
            P_array = np.concatenate((P_array, np.array([Pb])), axis=0)
            Y_array = np.concatenate((Y_array, np.array([Yb])), axis=0)
            Cap_array = np.concatenate((Cap_array, np.array([cap])))
            t_compute_array.append(t1-t0)

            print(f'b:{b:.3e}\tF:{Fb:.4f}')
            b = b * b_grow

        return F_array, b_array, P_array, Y_array, t_compute_array, Cap_array


    def anneal_DA_C(self, Y0, b_min, b_max, b_grow, b1_min, b1_max, b1_grow, theta, n_iters, D_tol):
       
        b = b_min
        Yb = Y0
        Db = self.get_D(Yb)
        Pb = self.get_PY_gibbs(Db, b)[0]

        F_array = []
        b_array = []
        P_array = np.array([Pb])
        Y_array = np.array([Yb])
        t_compute_array = []
        Cap_array = np.array([self.get_cap_cons(Pb, theta)[0] + self.C])

        while b <= b_max:
            t0 = time.time()
            if b < 1:
                b1_min = 0.01*b
                b1_max = b
                b1_grow = 1.2
            elif b > 1:
                b1_min = 0.01*b
                b1_max = b
                b1_grow = 1.2
            Pb, Yb = self.optimize_DA_C(Pb, Yb, b, b1_min, b1_max, b1_grow, theta, n_iters, D_tol)
            t1 = time.time()
            Fb = self.get_F(Yb, b)
            cap = self.get_cap_cons(Pb, theta)[0] + self.C
            Yb += np.random.multivariate_normal(np.zeros(2*self.M), 0.001*np.eye(2*self.M)).reshape(-1,2)
            b = b * b_grow

            F_array.append(Fb)
            b_array.append(b)
            P_array = np.concatenate((P_array, np.array([Pb])), axis=0)
            Y_array = np.concatenate((Y_array, np.array([Yb])), axis=0)
            Cap_array = np.concatenate((Cap_array, np.array([cap])))
            t_compute_array.append(t1-t0)

            print(f'b:{b:.3e}\tb1_min,b1_max:{b1_min, b1_max}\tF:{Fb:.4f}')
        
        return F_array, b_array, P_array, Y_array, t_compute_array, Cap_array
        

    ## ================= Clustering CBF ================= ##

    def h_hdot(self, P_ylx, u_p):
        h = P_ylx * (1 - P_ylx)
        hdot = cp.multiply((1 - 2 * P_ylx), u_p)
        return h, hdot

    def l_ldot(self, P_ylx, u_p):
        c = self.C
        q = self.rho.flatten()
        w = cp.multiply(P_ylx, self.alloc_C)
        w_dot = cp.multiply(u_p, self.alloc_C)
        l = c - q @ w
        ldot = - (q[:,np.newaxis].T @ w_dot).flatten()
        return l, ldot

    def s_sdot(self, P_ylx, u_p):
        b = self.B
        q = self.rho.flatten()
        w = cp.multiply(P_ylx, self.alloc_C)
        w_dot = cp.multiply(u_p, self.alloc_C)
        s = q @ w - b
        sdot = (q[:,np.newaxis].T @ w_dot).flatten()
        return s, sdot

    def F_gradF(self, x, beta):
        # compute cost matrix
        X = self.resLoc
        Y = x2Y(self.M, x)
        D = cdist(X, Y, metric = 'sqeuclidean')
        # free energy
        P = x2P(self.M,x)
        F = np.sum(self.rho*np.sum(P*(D + 1/beta * myLog(P)),axis=1,keepdims=True)) + 1/beta * np.log(self.M)
        # gradients
        P1 = np.expand_dims(P.T,axis=2)
        X_Y = np.array([y-self.resLoc for y in Y])
        dF_Y = 2*np.sum(self.rho * P1 * X_Y, axis=1)
        dF_P = self.rho*(D + 1/beta * (myLog(P) + 1))    
        dFdx = np.concatenate((dF_Y.flatten(),dF_P.flatten()))
        return F, dFdx

    def control_dyn(self, x, beta, u_b, p1, p2, gamma, alpha_h, alpha_l):
        N = self.N
        m = self.M
        q = self.rho.flatten()
        c = self.C
        Y = x2Y(self.M,x)
        P_ylx = x2P(self.M,x)
        # Decision variables 
        u_p = cp.Variable((N,m))
        u_y = cp.Variable((m,2))  
        delta = cp.Variable(1)
        # delta1 = cp.Variable(1)

        # Objective: minimize the sum of squares of x and the sum of q
        # objective = cp.Minimize(cp.sum_squares(u_p) + cp.sum_squares(u_y) + p1 * delta**2 + p2 * delta1**2)
        objective = cp.Minimize(cp.sum_squares(u_p) + cp.sum_squares(u_y) + p1 * delta**2)

        # Define constraints
        F, dF = self.F_gradF(x, beta)
        Fdot = cp.sum(cp.multiply(x2Y(self.M, dF), u_y)) + cp.sum(cp.multiply(x2P(self.M, dF), u_p))
        h, h_dot = self.h_hdot(P_ylx, u_p)
        l, l_dot = self.l_ldot(P_ylx, u_p)
        s, s_dot = self.s_sdot(P_ylx, u_p)

        constraints = [
            Fdot <= -gamma * F + delta,
            cp.sum(u_p, axis=1) == 0,
            h_dot >= -alpha_h * h,
            l_dot >= -alpha_l * l, #+ delta1
            s_dot >= -alpha_l * s #+ delta1
        ]

        # Define the problem
        problem = cp.Problem(objective, constraints)

        # Solver Options for OSQP
        solver_options = {
            'max_iter': 50000,         # Increase max iterations to 20000
            'eps_abs': 1e-4,           # Adjust absolute tolerance
            'eps_rel': 1e-4,           # Adjust relative tolerance
            'eps_prim_inf': 1e-3,      # Adjust primal infeasibility tolerance
            'eps_dual_inf': 1e-3,      # Adjust dual infeasibility tolerance
            'verbose': False           # Enable verbose output to track solver progress
        }

        # Solve the problem using OSQP with customized options
        result = problem.solve(solver = 'OSQP', **solver_options)

        # print(u_p.value)
        # print(u_y.value)
        
        # Check the results
        if np.isnan(problem.value).any() == True:
            print("Nan encountered!")
            return np.zeros((N,m)), np.zeros((m,2)), F, 0
        else:         
            return u_p.value, u_y.value, F, Fdot.value
    
    def dynamics(self, x, beta, p1, p2, gamma, alpha_h, alpha_l):    
        # Unpack the state vector
        Y = x2Y(self.M,x)
        P = x2P(self.M,x)
        u_b = 0
        # Compute control inputs; control_dyn should return u_p and u_y for state evolution.
        u_p, u_y, F_b, Fdot_b = self.control_dyn(x, beta, u_b, p1, p2, gamma, alpha_h, alpha_l) 
        # Flatten the derivatives into a single vector
        dxdt = YP2x(u_y, u_p)
        return dxdt, F_b, Fdot_b


    def optimize_CBF_CLF(self, x0, beta, p1, p2, gamma, alpha_h, alpha_l, T_f, dt_init=0.01, dt_min=1e-4, dt_max=0.1, Ftol=1e-2, xtol=1e-2, allowPrint=False):
        # initialization
        x_prev = x0
        dt_prev = dt_init 
        theta_prev = np.inf 
        t = 0.0
        iter_count = 0
        F_prev = 0.0

        while t < T_f:
            # Compute gradient step
            dxdt, F, Fdot = self.dynamics(x_prev, beta, p1, p2, gamma, alpha_h, alpha_l)

            # Compute new step size (lambda_k)
            if iter_count > 0:
                step_size_1 = np.sqrt(1 + theta_prev) * dt_prev
                grad_diff = np.linalg.norm(dxdt - dxdt_old) + 1e-6  # Regularization term to prevent division by zero
                step_size_2 = np.linalg.norm(x_prev - x_old) / (2 * grad_diff)  
                dt = min(max(step_size_2, dt_min), dt_max)  # Keep dt in range [dt_min, dt_max]
            else:
                dt = dt_init  # First iteration uses initial step size

            # Euler update
            x_next = x_prev + dt * dxdt

            # Projection step for P
            x_projected = YP2x(
                x2Y(self.M, x_next), 
                project_to_stochastic_matrix(x2P(self.M, x_next)))

            # Compute new theta_k
            if iter_count > 0:
                theta = dt / dt_prev
            else:
                theta = 1.0  # Initial theta value

            # Termination condition if change in F or x is too small
            if abs(F-F_prev) < Ftol:
                if allowPrint:
                    print(f"Ftol successful\titer:{iter_count}\ttime {t:.3e}\tFtol={abs(F-F_prev):.3e} < Ftol={Ftol:.3e}")
                break            
            elif allowPrint:
                print(f't{t:.3e}\tF:{F:.6f}\tFdiff:{abs(F-F_prev):6f}\tdt:{dt:.5f}\tdxdt_norm:{np.max(np.abs(dxdt)):.6f}\tdx:{np.max(np.abs(dxdt))*dt:.6f}\tsmallest_P:{np.min(x2P(self.M, x_prev))}')
            # if abs(Fdot*dt) < Ftol:
            #     if allowPrint:
            #         print(f"Ftol successful\titer:{iter_count}\ttime {t:.3e}\tFtol={abs(Fdot*dt):.3e} < Ftol={Ftol:.3e}")
            #     break
            # if max(abs(dxdt)*dt) < xtol:
            #     if allowPrint:
            #         print(f"xtol successful\titer:{iter_count}\ttime {t:.3e}\txtol={max(abs(dxdt*dt)):.3e} < xtol={xtol:.3e}")
            #     break
            
            # Update variables for next iteration
            x_old = x_prev
            dxdt_old = dxdt
            dt_prev = dt
            theta_prev = theta
            x_prev = x_projected
            t += dt
            iter_count += 1
            F_prev = F

        return x_prev, F

    # function to perform deterministic annealing
    def anneal_CBF_CLF(self, x0, beta0, betaf, beta_grow, p1, p2, gamma, alpha_h, alpha_l, T_f=15, dt_init=0.01, dt_min=1e-4, dt_max=0.1, Ftol=1e-2, xtol=1e-2, allowPrint=False, printOptimizeIters=False):
        xb = x0
        beta = beta0
        F_array = []
        b_array = []
        t_compute_array = []
        C_array = []
        P_array = np.array([x2P(self.M,x0)])
        Y_array = np.array([x2Y(self.M,x0)])
        Cap_array = np.array([np.sum(self.rho * self.alloc_C * x2P(self.M, x0),axis=0)])
        
        while beta <= betaf:
            
            # solve for the states
            t0 = time.time()
            xb, Fb = self.optimize_CBF_CLF(xb, beta, p1, p2, gamma, alpha_h, alpha_l, T_f, dt_init, dt_min, dt_max, Ftol, xtol, printOptimizeIters)
            t1 = time.time()
            # compute capacity
            cap = np.sum(self.rho * self.alloc_C * x2P(self.M, xb),axis=0)
            # perturb state
            xb[:2*self.M] += np.random.multivariate_normal(np.zeros(2*self.M), 0.001*np.eye(2*self.M))
            # update beta
            beta *= beta_grow
        
            # store data
            F_array.append(Fb)
            b_array.append(beta)
            P_array = np.concatenate((P_array, np.array([x2P(self.M,xb)])), axis=0)
            Y_array = np.concatenate((Y_array, np.array([x2Y(self.M,xb)])), axis=0)
            Cap_array = np.concatenate((Cap_array, np.array([cap])))
            t_compute_array.append(t1-t0)

            if allowPrint:
                print(f'beta:{beta:.3e}\tF:{Fb:.3e}')

        return F_array, b_array, P_array, Y_array, t_compute_array, Cap_array



    ## ================= Clustering OTS ================= ##

    # function to evaluate inequality constraints
    def G_ineq(self, x):
        P = x2P(self.M,x)
        H1 = ((P-self.P_eps)*(1-P)).flatten() # NM x 1
        H2 = self.C - np.sum(self.rho * self.alloc_C * P , axis=0)  # M x 1
        H3 = np.sum(self.rho * self.alloc_C * P ,axis=0) - self.B # M x 1
        G = np.concatenate((H1,H2,H3))
        return G

    # function to evaluate equality constraints
    def H_eq(self, x):
        P = x2P(self.M,x)
        H = np.sum(P,axis=1)-1
        return H

    # function to evaluate gradient of inequality constraints
    def gradG(self, x):
        N = self.N
        M = self.M
        P = x2P(self.M,x)
        dG11 = np.zeros((N*M, 2*M))
        dG12 = np.diag(1 + self.P_eps - 2*P.flatten())
        dG21 = np.zeros((M, 2*M))
        # dG22 = -np.kron((self.rho).T, np.eye(M))
        dG22 = -np.concatenate(list(map(np.diag, self.rho * self.alloc_C))).T
        dG31 = np.zeros((M, 2*M))
        # dG32 = np.kron((self.rho).T, np.eye(M))
        dG32  = np.concatenate(list(map(np.diag, self.rho * self.alloc_C))).T
        # print(f'shapes: dG11 {dG11.shape}, dG12 {dG12.shape}, dG21 {dG21.shape}, dG22 {dG22.shape}')
        dG = np.block([[dG11, dG12],[dG21, dG22], [dG31, dG32]])
        return dG

    # function to evaluate gradient of equality constraints
    def gradH(self, x):
        N = self.N
        M = self.M
        H1 = np.zeros((N,2*M))
        H2 = block_diag(*([np.ones((1,M))]*N))
        # print(f'shapes H1 {H1.shape} H2 {H2.shape}')
        dH = np.concatenate((H1,H2),axis=1)
        return dH


    # function to optimize free energy at given beta using SLSQP method
    def optimize_SLSQP(self, x0, beta, ftol):
        inequality_constraint = {'type':'ineq', 'fun':self.G_ineq, 'jac':self.gradG}
        equality_constraint = {'type':'eq', 'fun':self.H_eq, 'jac':self.gradH}
        F = lambda x, beta : self.F_gradF(x, beta)[0]
        gradF = lambda x, beta : self.F_gradF(x, beta)[1]
        res = minimize(
            F, x0,
            jac = gradF,
            args=(beta,), 
            method='SLSQP',
            constraints=[inequality_constraint, equality_constraint],
            options={'disp':False, 'ftol':ftol})
        return res.x, res.fun


    # function to perform deterministic annealing
    def anneal_SLSQP(self, x0, beta0, betaf, beta_grow, SLSQP_ftol):
        xb = x0
        beta = beta0
        F_array = []
        b_array = []
        P_array = np.array([x2P(self.M,x0)])
        Y_array = np.array([x2Y(self.M,x0)])
        t_compute_array = []
        Cap_array = np.array([np.sum(self.rho * self.alloc_C * x2P(self.M, x0),axis=0)])

        while beta <= betaf:
            t0 = time.time()
            xb, Fb = self.optimize_SLSQP(xb, beta, SLSQP_ftol)
            t1 = time.time()
            # compute capacity
            cap = np.sum(self.rho * self.alloc_C * x2P(self.M, xb),axis=0)
            # perturb states
            xb[:2*self.M] += np.random.multivariate_normal(np.zeros(2*self.M), 0.001*np.eye(2*self.M))
            # update beta
            beta *= beta_grow
        
            # store data
            F_array.append(Fb)
            b_array.append(beta)
            P_array = np.concatenate((P_array, np.array([x2P(self.M,xb)])), axis=0)
            Y_array = np.concatenate((Y_array, np.array([x2Y(self.M,xb)])), axis=0)
            Cap_array = np.concatenate((Cap_array, np.array([cap])))
            t_compute_array.append(t1-t0)

            # print
            print(f'beta:{beta:.3e}\tF:{Fb:.3e}')

        return F_array, b_array, P_array, Y_array, t_compute_array, Cap_array

    

    ## ================= clustering SGF ================= ##

    def calc_u1(self, x, beta, alpha, kappa):
        N = self.N
        M = self.M
        n_ineq = N*M + 2*M # number of inequality constraints
        n_eq = N # number of equality constraints

        # define variables and compute gradients
        u = cp.Variable(n_ineq + n_eq)
        dF = self.F_gradF(x, beta)[1]
        dG = self.gradG(x)
        dH = self.gradH(x)

        # cost and objective
        # objective = cp.Minimize(cp.sum_squares(u))
        objective = cp.Minimize(cp.sum_squares(dG.T @ u[0:n_ineq] + dH.T @ u[n_ineq:]))

        # constraints
        constraints = [
            (dG @ dG.T) @ u[0:n_ineq] + (dG @ dH.T) @ u[n_ineq:] <= -kappa * dG @ dF + alpha * self.G_ineq(x),
            u[0:N*M+M] <= np.zeros(N*M+M),
            -(dH @ dG.T) @ u[0:n_ineq] - (dH @ dH.T) @ u[n_ineq:] == kappa * dH @ dF - alpha * self.H_eq(x)
        ]

        # Solver Options
        solver_options = {
        'max_iter': 50000,         # Increase max iterations to 20000
        'eps_abs': 1e-4,           # Adjust absolute tolerance
        'eps_rel': 1e-4,           # Adjust relative tolerance
        'eps_prim_inf': 1e-3,      # Adjust primal infeasibility tolerance
        'eps_dual_inf': 1e-3,      # Adjust dual infeasibility tolerance
        'verbose': False            # Enable verbose output to track solver progress
        }

        # solve problem
        prob = cp.Problem(objective, constraints)
        prob.solve(solver='OSQP', **solver_options)

        return u.value[0:n_ineq], u.value[n_ineq:]
    
    def xdot(self, x, beta, alpha, kappa):
        # free energy
        F, dF = self.F_gradF(x, beta)
        # gradients
        # dF = self.gradF(x, beta)
        dG = self.gradG(x)
        dH = self.gradH(x)
        # control
        u, v = self.calc_u1(x, beta, alpha, kappa)
        # time derivatives
        dxdt = -kappa * dF - dG.T @ u - dH.T @ v
        dFdt = dF @ dxdt
        return dxdt, F, dFdt


    # function to solve ode at given beta
    def optimize_SGF(self, x0, beta, alpha, kappa, T_f, dt_init=0.01, dt_min=1e-4, dt_max=0.1, Ftol=1e-2, xtol=1e-2, allowPrint=False):
        N = self.N
        M = self.M
        x_prev = x0
        dt_prev = dt_init 
        theta_prev = np.inf 
        t = 0.0
        iter_count = 0
        F_prev = 0.0

        while t < T_f:
            dxdt, F, Fdot = self.xdot(x_prev, beta, alpha, kappa)

            # Compute new step size (lambda_k)
            if iter_count > 0:
                step_size_1 = np.sqrt(1 + theta_prev) * dt_prev
                grad_diff = np.linalg.norm(dxdt - dxdt_old) + 1e-6  # Regularization term to prevent division by zero
                step_size_2 = np.linalg.norm(x_prev - x_old) / (2 * grad_diff)  
                dt = min(max(step_size_2, dt_min), dt_max)  # Keep dt in range [dt_min, dt_max]
            else:
                dt = dt_init  # First iteration uses initial step size

            # Euler update
            x_next = x_prev + dt * dxdt

            # Projection step for P
            x_projected = YP2x(
                x2Y(self.M, x_next), 
                project_to_stochastic_matrix(x2P(self.M, x_next)))

            # Compute new theta_k
            if iter_count > 0:
                theta = dt / dt_prev
            else:
                theta = 1.0  # Initial theta value

            # Termination condition if |Fdot| is too small
            if abs(F-F_prev) < Ftol:
                if allowPrint:
                    print(f"Ftol successful\titer:{iter_count}\ttime {t:.3e}\tFtol={abs(F-F_prev):.3e} < Ftol={Ftol:.3e}")
                break            
            # if abs(Fdot*dt) < Ftol:
            #     if allowPrint:
            #         print(f"Ftol successful\titer:{iter_count}\ttime {t:.3e}\tFtol={abs(Fdot*dt):.3e} < Ftol={Ftol:.3e}")
            #     break
            # if max(abs(dxdt)*dt) < xtol:
            #     print(f"xtol successful\titer:{iter_count}\ttime {t:.3e}\txtol={max(abs(dxdt*dt)):.3e} < xtol={xtol:.3e}")
            #     break
            elif allowPrint:
                print(f't{t:.3e}\tF:{F:.6f}\tFdiff:{abs(F-F_prev):6f}\tdF:{Fdot*dt:.3f}\tdt:{dt:.5f}\tdxdt_norm:{np.max(np.abs(dxdt)):.6f}\tdx:{np.max(np.abs(dxdt))*dt:.6f}')
            
            # Update variables for next iteration
            x_old = x_prev
            dxdt_old = dxdt
            dt_prev = dt
            theta_prev = theta
            x_prev = x_projected
            t += dt
            iter_count += 1
            F_prev = F
            
        return x_prev, F

    # function to perform deterministic annealing
    def anneal_SGF(self, x0, beta0, betaf, beta_grow, alpha, kappa, T_f=15, dt_init=0.01, dt_min=1e-4, dt_max=0.1, Ftol=1e-2, xtol=1e-2, allowPrint=False, printOptimizeIters=False):
        xb = x0
        beta = beta0
        F_array = []
        b_array = []
        P_array = np.array([x2P(self.M, x0)])
        Y_array = np.array([x2Y(self.M, x0)])
        t_compute_array = []
        Cap_array = np.array([np.sum(self.rho * self.alloc_C * x2P(self.M, x0),axis=0)])

        while beta <= betaf:
            t0 = time.time()
            xb, Fb = self.optimize_SGF(xb, beta, alpha, kappa, T_f, dt_init, dt_min, dt_max, Ftol, xtol, allowPrint=printOptimizeIters)
            t1 = time.time()
            # compute capacity
            cap = np.sum(self.rho * self.alloc_C * x2P(self.M, xb),axis=0)
            # perturb state
            xb[:2*self.M] += np.random.multivariate_normal(np.zeros(2*self.M), 0.001*np.eye(2*self.M))
            # update beta
            beta *= beta_grow
        
            # store data
            F_array.append(Fb)
            b_array.append(beta)
            P_array = np.concatenate((P_array, np.array([x2P(self.M, xb)])), axis=0)
            Y_array = np.concatenate((Y_array, np.array([x2Y(self.M, xb)])), axis=0)
            Cap_array = np.concatenate((Cap_array, np.array([cap])))
            t_compute_array.append(t1-t0)

            if allowPrint:
                print(f'beta:{beta:.3e}\tF:{Fb:.3e}')

        return F_array, b_array, P_array, Y_array, t_compute_array, Cap_array
        
    
# function to plot FLP results
def plot_flp(flp, res_means, Y_arr, P_arr, figSize):
    ''' 
    input - flp : the class instance
            N : # resources
            M : # facilities
            resLoc : resource locations
            facLoc0 : initial facility locations
            facLoc_opt : final facility locations (supposedly optimal)
            ineq_figures : data for plotting inequality constraints, datatype is variable
    output: None
    '''
    N = flp.N
    M = flp.M
    resLoc = flp.resLoc
    Y_init = Y_arr[0]
    P_init = P_arr[0]
    Y_final = Y_arr[-1]
    P_final = P_arr[-1]
    # resource and data colours
    np.random.seed(3)
    resource_colors = np.random.uniform(0,1,(flp.M, 3))
    resource_colors = resource_colors/resource_colors.sum(axis=1,keepdims=True)
    data_colors = np.dot(P_final, resource_colors)
    # initialize figure
    plt.figure(figsize=figSize)
    # demand point locations
    plt.scatter(resLoc[:, 0], resLoc[:, 1], marker="o", c=data_colors, edgecolor='k', s=100, alpha=0.8, label="Demand Points")
    # demand point means
    plt.scatter(res_means[:, 0], res_means[:, 1], marker = ".")
    # initial resource locations
    plt.scatter(Y_init[:, 0], Y_init[:, 1], marker = "x", c=resource_colors, s=100, edgecolor='k', alpha=0.4, label='Initial Resource Locations')
    # final resource locations
    plt.scatter(Y_final[:, 0], Y_final[:, 1], c=resource_colors, marker = "*", s=200, edgecolor='yellow', label='Final Resource Positions')
    # centroid locations
    if len(Y_arr) > 1:
        Y_centroid = Y_arr[1]
        plt.scatter(Y_centroid[:, 0], Y_centroid[:, 1], c=resource_colors, marker = "s", s=100, label=rf'Resource Positions at $\beta=\beta_{{min}}$')
        
    # Plot trajectories and final positions of resources
    for i in range(M):
        plt.plot(
            Y_arr[:, i, 0], Y_arr[:, i, 1], linestyle="dotted", linewidth=2, color=resource_colors[i])
        if i == M-1:
            plt.plot(
                Y_arr[:, i, 0], Y_arr[:, i, 1],
                label=f"Resource Trajectories", linestyle="dotted", linewidth=2, color=resource_colors[i])
        
    # final facility locations
    fontSize = 18
    # plt.title(f'(N,M)={N,M}', fontsize=fontSize)
    plt.legend(fontsize=0.7*fontSize)
    plt.xticks(fontsize=fontSize)
    plt.yticks(fontsize=fontSize)
    plt.grid()
    plt.show()
    return None


def plot_flp_with_capacity_1(flp, res_means, Y_arr, P_arr, 
                x, data, min_vals, max_vals, figSize, fontSize, 
                x_expand_factor, y_expand_factor, filedir, 
                isSavePDF=False, inset_size=(1, 3), 
                inset_loc=[0.79, 0.2, 0.1, 0.35], plotResTraj=False):
    ''' 
    input - flp : class instance containing resource and facility data
            res_means, Y_arr, P_arr : arrays for resource means and trajectories
            x, data, min_vals, max_vals : data for the inset plot
            figSize : main figure size
            fontSize : font size for labels
            inset_size : tuple defining the inset size (width, height) as a fraction of main plot
            inset_loc : location of the inset ('upper right', 'lower left', etc.)
    output: None
    '''
    
    # Extract common data
    N = flp.N
    M = flp.M
    resLoc = flp.resLoc

    # Generate colors for resources
    np.random.seed(4)
    resource_colors = np.random.uniform(0,1,(M, 3))
    resource_colors = resource_colors/resource_colors.sum(axis=1, keepdims=True)

    # Extract solution specific data
    Y_init = Y_arr[0]
    P_final = P_arr[-1]
    Y_final = Y_arr[-1]

    # solution specific demand point colors
    data_colors = np.dot(P_final, resource_colors)

    # Create the main figure
    fig, ax = plt.subplots(figsize=figSize)
    
    # Main scatter plot (facility locations)
    ax.scatter(resLoc[:, 0], resLoc[:, 1], marker="X", c=data_colors, edgecolor='white', s=600, alpha=0.3, label="Users")
    ax.scatter(res_means[:, 0], res_means[:, 1], marker = ".")

    # Scatter for Y_final with labels
    ax.scatter(Y_final[:, 0], Y_final[:, 1], c=resource_colors, marker="s", s=1000, edgecolor='yellow', label='R final')
    for i in range(M):
        ax.text(Y_final[i, 0], Y_final[i, 1], str(i), fontsize=fontSize*0.8, color='yellow', ha='center', va='center')

    # Plot resource trajectories
    if plotResTraj:
        # Scatter for Y_init with labels
        ax.scatter(Y_init[:, 0], Y_init[:, 1], marker=".", c=resource_colors, s=600, alpha=0.8, label='R initial')
        # for i in range(M):
        #     ax.text(Y_init[i, 0] + 0.2, Y_init[i, 1] + 0.2, str(i), fontsize=fontSize*0.8, color='black')
        # centroid locations
        if len(Y_arr) > 1:
            Y_centroid = Y_arr[1]
            plt.scatter(Y_centroid[:, 0], Y_centroid[:, 1], c=resource_colors, marker = "s", s=100, label=rf'R at $\beta={0.001}$')

        for i in range(M):
            ax.plot(Y_arr[:, i, 0], Y_arr[:, i, 1], linestyle="dotted", linewidth=2, color=resource_colors[i])
            if i == M-1:
                ax.plot(Y_arr[:, i, 0], Y_arr[:, i, 1], label="R Trajectory", linestyle="dotted", linewidth=4, color=resource_colors[i])

    # Expand x-axis limits to create room for the inset
    xmin, xmax = ax.get_xlim()
    x_range = xmax - xmin
    new_xmax = xmax + x_expand_factor * x_range  # Expand xmax by a fraction of x-range
    ax.set_xlim(xmin, new_xmax)

    # # Expand x-axis limits to create room for the inset
    # ymin, ymax = ax.get_ylim()
    # y_range = ymax - ymin
    # newca_ymax = ymax + y_expand_factor * y_range  # Expand xmax by a fraction of x-range
    # ax.set_xlim(ymin, new_ymax)

    # Format main plot
    ax.legend(fontsize=0.7*fontSize, ncol=1, loc="upper right")
    ax.tick_params(axis='both', labelsize=fontSize)
    ax.grid()

    # Add an inset axes for the allocation cost plot (placed at bottom right)
    inset_ax = fig.add_axes(inset_loc)
    # inset_ax = inset_axes(ax, width=inset_size[0], height=inset_size[1], loc=inset_loc)
    
    # "Transpose" the inset plot: Make the x-axis "Allocation Cost" and y-axis "Resource Index"
    
    # Plot bars for the min-max range (horizontal instead of vertical)
    for i in range(len(x)):
        inset_ax.barh(x[i], max_vals[i] - min_vals[i], left=min_vals[i], color='lightgray', height=0.2, alpha=1.0)
    
    # Scatter plot for actual data points (horizontal instead of vertical)
    inset_ax.scatter(data, x, color='black', marker='o', s=50, label='Final Allocation Cost')
    
    # Format inset plot
    inset_ax.set_ylabel("R ID", fontsize=fontSize*0.8)  # Now y-axis represents resource index
    inset_ax.set_yticks([i for i in range(flp.M)])
    inset_ax.set_xlabel("U", fontsize=fontSize*0.8)  # Now x-axis represents allocation cost
    inset_ax.tick_params(axis='both', labelsize=fontSize*0.8)
    # inset_ax.legend(fontsize=fontSize*0.5)


    if isSavePDF:
        fig.savefig(filedir, format="pdf")
    
    plt.show()


# def plot_data_with_shaded_region(x, data, min_vals, max_vals, figSize, fontSize):
#     plt.figure(figsize=figSize)
    
#     # Plot bars for the min-max range at each index
#     for i in range(len(x)):
#         plt.bar(x[i], max_vals[i] - min_vals[i], bottom=min_vals[i], color='lightgray', width=0.5, alpha=0.8)
    
#     # Scatter plot for actual data points
#     plt.scatter(x, data, color='red', marker='_', s=1000, label='Final Allocation Cost')
    
#     plt.xlabel("Resource Index", fontsize=fontSize*0.8)
#     plt.ylabel("Allocation Cost", fontsize=fontSize*0.8)
#     plt.xticks(fontsize=fontSize)
#     plt.yticks(fontsize=fontSize)
#     # plt.title("Data with Min-Max Range", fontsize=fontSize)
#     plt.legend(fontsize=fontSize)
#     plt.show()


