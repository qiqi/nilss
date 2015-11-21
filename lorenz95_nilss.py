from numpy import *
from lorenz95 import dot, step, adjoint, gradContribution
from nilss import NILSS

s = 20
nChunks, nStepsPerChunk = 200, 50
nTotalSteps = nChunks * nStepsPerChunk

# initial transient
x = random.rand(40)
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
nHomo = 20
lss = NILSS(nHomo, x.size, dot)

y = random.rand(nHomo + 1, x.size)
for iChunk in reversed(range(nChunks)):
    grad = zeros(y.shape[0])
    for iStep in reversed(range(nStepsPerChunk)):
        x = history[iChunk][iStep]
        for iAdj in range(y.shape[0]):
            grad[iAdj] += gradContribution(x, y[iAdj])
            y[iAdj] = adjoint(x, y[iAdj], s)
        timeFraction = (iStep + iChunk * nStepsPerChunk) / nTotalSteps
        windowFunction = sin(timeFraction * pi)**2 * 2
        y[-1] += windowFunction / nTotalSteps # inhomogeneous
    lss.checkpoint(y, grad)
    print(lss.grad())
