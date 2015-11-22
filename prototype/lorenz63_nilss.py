from numpy import *
from lorenz63 import dot, step, adjoint, gradContribution
from nilss import NILSS

s = 28
nChunks, nStepsPerChunk = 100, 1000
nTotalSteps = nChunks * nStepsPerChunk

# initial transient
x = random.rand(3)
for iStep in range(1000):
    x = step(x, s)

# primal
history = []
for iChunk in range(nChunks):
    history.append([])
    for iStep in range(nStepsPerChunk):
        history[iChunk].append(x)
        x = step(x, s)

# adjoint
nHomo = 2
lss = NILSS(nHomo, x.size, dot)

y = random.rand(nHomo + 1, x.size)
for iChunk in reversed(range(nChunks)):
    grad = zeros(y.shape[0])
    yHist = []
    for iStep in reversed(range(nStepsPerChunk)):
        x = history[iChunk][iStep]
        for iAdj in range(y.shape[0]):
            grad[iAdj] += gradContribution(x, y[iAdj])
            y[iAdj] = adjoint(x, y[iAdj], s)
        timeFraction = (iStep + iChunk * nStepsPerChunk) / nTotalSteps
        windowFunction = sin(timeFraction * pi)**2 * 2
        y[-1][2] += windowFunction / nTotalSteps # inhomogeneous
        yHist.append(y.copy())
    lss.checkpoint(y, grad, yHist)
    print(lss.grad())

y = array(lss.y_hist)
a = array(lss.a)[:-1]
y = (y[:,:,:-1,:] * a[:,newaxis,:,newaxis]).sum(2) + y[:,:,-1,:]
y = y.reshape((-1,3))[::-1,:]

x = array(history).reshape([-1,3])

g = x[:,0] * y[:,1] * 0.001