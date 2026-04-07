import numpy as np
import community_sim

def one_comm_dynamics(t, N, a, b, alpha):
    dNdt = N*(b*(1 - alpha*N) + (a@N))
    return dNdt

def one_comm_dynamics2(t, N, a0, b):
    dNdt = N*(b + (a0@N))
    return dNdt

t = 0
N = np.array([10,20,30])
alpha = 0.05
b = np.array([0.01, 0.01, 0.01])

a = np.array([[0, 0.005, 0.005],
              [0.005, 0, 0.005],
              [0.005, 0.005, 0]])

a0 = np.array([[-alpha * b[0], 0.005, 0.005],
              [0.005, -alpha * b[1], 0.005],
              [0.005, 0.005, -alpha * b[2]]])

dNdt0 = community_sim.one_comm_dynamics(t, N, a, b, alpha)
dNdt1 = one_comm_dynamics(t, N, a, b, alpha)
dNdt2 = one_comm_dynamics2(t, N, a0, b)

def test_one_comm_dynamics():
    assert(all(dNdt0 == dNdt1))
    assert(all(dNdt1 == dNdt2))

