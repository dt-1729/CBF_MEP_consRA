import numpy as np
from scipy.spatial.distance import cdist


def myLog(x):
    logx = np.log(np.maximum(1e-20,x))
    # logx = np.log(x)
    if np.isnan(logx).any() == True:
        print(f'Inside myLog: NaN encountered!\nmin_x:{x.min()}')
        return None
    else:
        return logx

def x2Y(M, x):
    Y = x[0:2*M].reshape(-1,2)
    return Y

def x2P(M, x):
    P = x[2*M:].reshape(-1,M)
    return P

def YP2x(Y, P):
    y = Y.flatten()
    p = P.flatten()
    x = np.concatenate((y,p))
    return x

# Auxiliary functions to project back into the feasible space
def project_to_stochastic_matrix(matrix):
    """
    Projects each row of the input matrix onto the simplex (triangle) defined by:
    - The row sums to 1
    - Each element in the row lies within [0, 1]

    :param matrix: np.ndarray, the input matrix (shape: NxM)
    :return: np.ndarray, the projected matrix (same shape as input)
    """
    def project_to_plane(vector):
        N = vector.shape[0]
        normal = np.ones((N,))
        p = 1 / N * np.ones((N,))
        n_dot_n = N
        v_dot_n = np.dot(vector, normal)
        p_dot_n = np.dot(p, normal)
        projection = vector - ((v_dot_n - p_dot_n) / n_dot_n) * normal
        return projection

    def project_to_triangle(v):
        tol = 1e-7
        v_proj = project_to_plane(v)
        v_clamped = np.clip(v_proj, tol, 1 - tol)
        sum_clamped = np.sum(v_clamped)

        if sum_clamped == 1:
            return v_clamped
        elif sum_clamped < 1:
            deficit = 1 - sum_clamped
            free_indices = v_clamped < 1
            num_free = np.sum(free_indices)
            if num_free > 0:
                increment = deficit / num_free
                v_clamped[free_indices] += increment
            return v_clamped
        else:
            return v_clamped / sum_clamped

    # Apply the projection to each row of the matrix
    projected_matrix = np.apply_along_axis(project_to_triangle, axis=1, arr=matrix)
    return projected_matrix

def generate_testcase(
    N, 
    M, 
    d, 
    cluster_cov, 
    seed=0, 
    sq_size=20, 
    unif_split=True, 
    normalizer=False):
    # seed for reproducibility
    np.random.seed(seed)
    # cluster means and covariances
    mean = [np.random.uniform(-sq_size,sq_size,d) for i in range(M)]
    cov = [cluster_cov*np.diag(np.random.uniform(0,1,(2,))) for i in range(M)]
    # cluster split percentage
    if unif_split:
        split_pct = np.ones(M)
    else:
        split_pct = np.random.uniform(0.01, 2, M)    
    split_pct = split_pct/split_pct.sum()

    resLoc = np.zeros([N, d]) # initialize resource locations
    res_mean = np.zeros([len(mean), d]) # initialize resource mean locations
    rho = np.ones([N,1])/N # resource weights/ uniform for now
    # generate resource and facility locations using normal distribution
    cnt = 0
    for i in range(len(mean)):
        n_split = int(split_pct[i]*N)
        resLoc[cnt:cnt+n_split, :] = np.random.multivariate_normal(mean[i], cov[i], n_split)
        res_mean[i,:] = np.sum(resLoc[cnt:cnt+n_split, :], axis=0)/len(resLoc[cnt:cnt+n_split, :])
        cnt = cnt + n_split
    facLoc = np.tile(np.sum(resLoc*rho, axis=0), (M,1)) # initialize facility locations at the centroid
    # compute normalizer
    resLoc_distances = cdist(resLoc, resLoc, metric='euclidean')
    if normalizer:
        normalizer = np.max(resLoc_distances)
    else:
        normalizer = 1.0

    return resLoc/normalizer, facLoc/normalizer, res_mean/normalizer, split_pct, rho

# def generate_multiple_testcases( 
#     N_range, 
#     M_range, 
#     cluster_cov_range, 
#     sq_size, 
#     d=2, 
#     normalizer=False):

#     for N, M, cov in zip(N_range, M_range, cluster_cov_range):
#         testcase = {}
#         X, Y0, Xm, split_pct, rho = generate_testcase(
#             N, 
#             M, 
#             d, 
#             cluster_cov, 
#             seed=0, 
#             sq_size=sq_size, 
#             unif_split=unif_split,
#             normalizer=normalizer)

#         alloc_cost = np.random.uniform(1,1,(N,M))
#         mean_alloc_cost = np.mean(alloc_cost)
#         C = np.random.uniform(0.3,1,len(split_pct))
#         # C = np.ones(shape=split_pct.shape)
#         C = (C/np.sum(C) + 0.02) * mean_alloc_cost
#         B = np.random.uniform(C.min()*0.3, C.min()*0.6, len(split_pct))
        
#         testcase = {
#             'X':X, 
#             'Y0':Y0,
#             'Xm':Xm,
#             'split_pct':split_pct,
#             'rho':rho
#         }
