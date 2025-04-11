import numpy as np
from scipy.spatial.distance import cdist



def generate_testcase(N, M, d, cluster_cov, seed=0, sq_size=20, unif_split=True, normalizer=False):
    # seed for reproducibility
    np.random.seed(seed)
    # cluster means and covariances
    mean = [np.random.uniform(-sq_size,sq_size,d) for i in range(M)]
    cov = [cluster_cov*np.diag(np.random.uniform(0,1,(2,))) for i in range(M)]
    # cluster split percentage
    if unif_split:
        split_pct = np.ones(M)
    else:
        split_pct = np.random.uniform(0.01,2,M)    
        # split_pct = np.array([0.2, 0.3, 0.5])
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



# function to initialize various testcases
def testcase_catalog(tc_name):
    ''' 
    input - tc_name: name of the testcase
    output - N: number of resources
             M: number of facilities
             resLoc: resource locations
             facLoc: facility locations (initial)
    '''
    if 'small_tc_160_4' == tc_name :
        d = 2 # dimension of the coordinates
        N = 160 # number of resources
        M = 4 # number of facilties
        mean = [np.array([-5,0]), np.array([15,20]), np.array([-20,-15]), np.array([-25,25])] # cluster means
        cov = [np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3] # cluster covariances
        split_pct = np.array([0.26, 0.20, 0.24, 0.30]) # no of points in every cluster
        C = np.array([0.4, 0.4, 0.4, 0.4]) # the capacity of the clusters; C <= split
    elif 'small_tc_320_8' == tc_name :
        d = 2 # dimension of the coordinates
        N = 320 # number of resources
        M = 8 # number of facilties
        mean = [np.array([-5,0]), np.array([15,20]), np.array([-20,-15]), np.array([-25,25]),
                np.array([15,0]), np.array([-15,-30]), np.array([40,15]), np.array([25,-25])] # cluster means
        cov = [np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3] # cluster covariances
        split_pct = np.random.uniform(5,10,(8,)) # no of points in every cluster
        split_pct = split_pct/np.sum(split_pct)
    elif 'small_tc_44_2' == tc_name :
        d = 2 # dimension of the coordinates
        N = 44 # number of resources
        M = 2 # number of facilties
        mean = [np.array([-15,0]), np.array([15,0])] # cluster means
        cov = [np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3] # cluster covariances
        split_pct = np.array([0.55, 0.45]) # no of points in every cluster
    elif 'small_tc_48_3' == tc_name :
        d = 2 # dimension of the coordinates
        N = 48 # number of resources
        M = 3 # number of facilties
        mean = [np.array([-10,0]), np.array([10,0]), np.array([0,10])] # cluster means
        cov = [np.array([[1, 0],[0, 1]])*3,
               np.array([[1, 0],[0, 1]])*3,
               np.array([[1, 0],[0, 1]])*3] # cluster covariances
        split_pct = np.array([0.3, 0.4, 0.3])
    elif 'small_tc_10_2' == tc_name:
        d = 2 # dimension of the coordinates
        N = 10 # number of resources
        M = 2 # number of facilties
        mean = [np.array([-15,0]), np.array([15,0])] # cluster means
        cov = [np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3] # cluster covariances
        split_pct = np.array([0.6, 0.4]) # no of points in every cluster
    elif 'small_tc_4_2' == tc_name:
        d = 2 # dimension of the coordinates
        N = 4 # number of resources
        M = 2 # number of facilties
        mean = [np.array([-5,0]), np.array([5,0])] # cluster means
        mean = mean
        cov = [np.array([[5, 0],[0, 5]])*3,
               np.array([[5, 0],[0, 5]])*3] # cluster covariances
        split_pct = np.array([0.5, 0.5]) # no of points in every cluster
    elif 'small_tc_300_3' == tc_name:
        d = 2 # dimension of the coordinates
        N = 300 # number of resources
        M = 3 # number of facilties
        mean = [np.array([-10,1]), np.array([-5,15]), np.array([9,10])] # cluster means
        mean = mean
        cov = [np.array([[1, 0],[0, 1]])*5,
               np.array([[1, 0],[0, 1]])*3,
               np.array([[1, 0],[0, 1]])*4] # cluster covariances
        split_pct = np.array([0.33, 0.33, 0.34]) # no of points in every cluster
    elif 'large_tc_1000_10' == tc_name:
        np.random.seed(2)
        d = 2 # dimension of the coordinates
        N = 1000 # number of resources
        M = 10 # number of facilties
        mean = [np.random.uniform(-20,20,2) for i in range(M)]
        # mean = [np.array([-10,1]), np.array([-5,15]), np.array([9,10])] # cluster means
        mean = mean
        cov = [np.eye(2)*np.random.uniform(0,5) for i in range(M)]
        split_pct = np.ones(M)
        # split_pct = np.random.uniform(0,1,M)
        split_pct = split_pct/split_pct.sum()

    # C = np.array([1, 1]) # the capacity of the clusters; C <= split
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
    normalizer = 1.0 #np.max(resLoc_distances)/10
    # print(res_mean)
    return N, M, d, resLoc/normalizer, facLoc/normalizer, res_mean/normalizer, split_pct, rho

