from numpy import *

dt = 0.001

def step(x, s):
    return dt * (-roll(x, 2) * roll(x, 1) + roll(x, 1) * roll(x, -1) + s) \
         + (1 - dt) * x

def adjoint(x, y, s):
    return dt * (-roll(x, 1) * roll(y, -1) - roll(x, -1) * roll(y, -2) \
               + roll(x, 2) * roll(y, 1) + roll(x, -2) * roll(y, -1)) \
         + (1 - dt) * y

def gradContribution(x, y):
    return y.sum() * dt

if __name__ == '__main__':
    nTests = 10
    for i in range(nTests):
        s = random.rand() * 100
        x = random.rand(40)
        y = random.rand(40)
        eps = 1E-8
        dx = random.rand(40) * eps/2
        grad1 = dot(step(x+dx, s) - step(x-dx, s), y) / eps
        grad2 = dot(2*dx, adjoint(x, y, s)) / eps
        print('FD = {0}, ADJ = {1}, error = {2}'.format(grad1, grad2,
            grad1 - grad2))

    for i in range(nTests):
        s = random.rand() * 100
        x0 = random.rand(40)
        y0 = random.rand(40)
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
        print('FD = {0}, ADJ = {1}, error = {2}'.format(grad1, grad2,
            grad1 - grad2))
