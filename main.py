import itertools
import argparse
import pyqubo
from pyqubo import Spin
import numpy
import sympy as sym
import math
import multiprocessing
from pulp import *

def sat1(x):
    return not x[0] and not x[1] and not x[2] and not x[3] and not x[4] and x[5] and not x[6] and x[7] and x[8] and not x[9] and x[10] and x[11] and x[12] and x[13] and x[14] and not x[15]

def convert_to_minus(x):
    # 0 -> 1, 1 -> -1
    return (1 - 2*numpy.multiply(x, 1))

N=16
J = numpy.array(sym.Matrix(N, N, lambda i,j: sym.Symbol('J_%d,%d' % (i, j))).tolist())

def getEquations(*x):
  a = numpy.zeros(N*(N+1)+1)
  for i in range(N):
    a[i] = convert_to_minus(x[i]) # find coefficients for h_i
  for row in range(len(J)):
    for col in range(len(J[row])):
      if (row < col):
        a[col + row*N + N] = convert_to_minus(x[row])*convert_to_minus(x[col]) # find coefficients for J_ij, i<j
      else:
        a[col + row*N + N] = 0
  sat = sat1(x)
  if (sat == True):
    a[N*(N+1)] = 0
    return a
  else:
    a[N*(N+1)] = -1 # find coefficients for g
    return a

def main():
  num_cores = multiprocessing.cpu_count()
  pool = multiprocessing.Pool(num_cores)
  res = pool.starmap_async(getEquations, list(itertools.product([False, True], repeat=16)), chunksize=10)
  pool.close()
  pool.join()
  eq_arr = res.get()

  prob = LpProblem("Hamiltonian_coeff", LpMinimize)
  h = LpVariable.dicts("h", range(0, N), lowBound=-2, upBound=2)
  JJ = LpVariable.dicts("JJ", range(N, N**2+N), lowBound=-1, upBound=1)
  g = LpVariable.dicts("g", range(N+N**2, N+N**2+1), lowBound=0)
  prob += (-1)*g[N+N**2]
  total_dict = {**h, **JJ, **g} # variables
  var_list = list(total_dict.values())
  print("var_list", var_list)
  prod = numpy.matmul(eq_arr,var_list)
  print("prod", prod)
  print("constraints")
  for i in range(len(prod)):
    if (eq_arr[i][(N+N**2)] == -1):
      print(LpConstraint(prod[i],rhs=0,sense = 1))
      prob += LpConstraint(prod[i],rhs=0,sense = 1)
    else:
      print(LpConstraint(prod[i],rhs=0,sense = 0))
      prob += LpConstraint(prod[i],rhs=0,sense = 0)
  status = prob.solve()

  print("status\n") 
  print(LpStatus[status])
  print("\nh\n")
  for i in (range(len(h))):
    print(h[i].varValue)
  print("\nJJ\n")
  for i in range(N):
    for j in range(N):
      print(JJ[i*N+j+N].varValue)
    print("\n")
  print("\ng\n")
  print(g[N+N**2].varValue)

  return 0



if __name__ == '__main__':
    main()
