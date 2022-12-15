import scipy.linalg as LA
from pyargus import directionEstimation as de
import numpy as np
b = np.array([[1,2,3,4,5],
              [1.1,2,3.5,4,6],
              [2,3,4,5,1.1],
              [2,3.5,1,4,6],
              [3,4,5.5,6,1],
              [3.1,4,4,5,6]])


CovMat = de.corr_matrix_estimate(b.T, imp="fast")
print(CovMat)

invR=LA.inv(CovMat)

print(invR*2)