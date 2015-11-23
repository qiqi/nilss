from __future__ import division
from numpy import *

dt = 0.001

def step(x, s):
    return x + dt * array([10 * (x[1] - x[0]),
                           x[0] * (s - x[2]) - x[1],
                           x[0] * x[1] - 8./3 * x[2]])

def adjoint(x, y, s):
    return y + dt * array([-10 * y[0] + (s - x[2]) * y[1] + x[1] * y[2],
                           10 * y[0] - y[1] + x[0] * y[2],
                           -x[0] * y[1] - 8./3 * y[2]])

def gradContribution(x, y):
    return x[0] * y[1] * dt

if __name__ == '__main__':
    nTests = 10
    for i in range(nTests):
        s = random.rand() * 100
        x = random.rand(3) * 50
        y = random.rand(3)
        eps = 1E-8
        dx = random.rand(3) * eps/2
        grad1 = dot(step(x+dx, s) - step(x-dx, s), y) / eps
        grad2 = dot(2*dx, adjoint(x, y, s)) / eps
        print('FD = {0}, ADJ = {1}, error = {2}'.format(grad1, grad2,
            grad1 - grad2))

    for i in range(nTests):
        s = random.rand() * 100
        x0 = random.rand(3) * 50
        y0 = random.rand(3)
        eps = 1E-8
        ds = eps / 2
        # primal finite difference
        xp, xm, x = [x0], [x0], [x0]
        for i in range(10):
            xp.append(step(xp[-1], s+ds))
            xm.append(step(xm[-1], s-ds))
            x.append(step(x[-1], s-ds))
        grad1 = (dot(xp[-1], y) - dot(xm[-1], y)) / eps
        # adjoint
        grad2 = 0
        for i in reversed(range(10)):
            grad2 += gradContribution(x[i], y)
            y = adjoint(x[i], y, s)
        print('FD={0}, ADJ={1}, error={2}'.format(grad1, grad2,
            grad1 - grad2))
