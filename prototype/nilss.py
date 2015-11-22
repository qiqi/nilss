from copy import copy as pythonCopy
from numpy import *
from scipy import sparse
from scipy.sparse import linalg as splinalg

_DEBUG_ = True
_DEBUG_ = False

def LSS_KKT(R, D):
    R, D = array(R), array(D)
    assert R.ndim == 3
    assert R.shape[1] == R.shape[2]
    N, m = R.shape[:2]

    bigR = sparse.bsr_matrix((R, r_[:N], r_[:N+1]), \
                          shape=(N*m, (N+1)*m))

    I = array([eye(m)] * N)
    bigI = sparse.bsr_matrix((I, r_[1:N+1], r_[:N+1]), \
                          shape=(N*m, (N+1)*m))

    bigL = bigI - bigR

    assert D.shape == (N+1, m, m)
    bigD = sparse.bsr_matrix((D, r_[:N+1], r_[:N+2]), \
                          shape=((N+1)*m, (N+1)*m))

    O = zeros([N, m, m])
    bigO = sparse.bsr_matrix((O, r_[:N], r_[:N+1]), \
                          shape=(N*m, N*m))

    return sparse.vstack([sparse.hstack([bigD, bigL.T]),
                          sparse.hstack([bigL, bigO])])

def solveLss(R, D, b, c):
    N, m = len(R), len(R[0])

    kkt = LSS_KKT(R, D)

    assert(len(R) == len(b))
    assert(len(D) == len(c))

    rhs = hstack([hstack(c), hstack(b)])

    sol = splinalg.spsolve(kkt, rhs)
    return sol[:(N+1)*m].reshape([N+1,m])


class NILSS(object):
    def __init__(self, nHomo, size, dot):
        self.nHomo = nHomo
        self.size = size
        self.dot = dot

        self.R = []
        self.b = []
        self.stored_grad = []

        self.y_last = []
        self.y_next = []

        self.computed_grad = []
        self.y_hist = []

    def checkpoint(self, y, grad, yHist=None):
        if _DEBUG_:
            self.y_last.append([yi.copy() for yi in y])
            if yHist:
                self.y_hist.append(yHist)

        cov = array([[self.dot(yi, yj) for yi in y[:-1]] for yj in y[:-1]])
        assert(cov.shape == (self.nHomo, self.nHomo))
        R = linalg.cholesky(cov).T
        self.R.append(R)

        for i in range(self.nHomo):
            for j in range(i):
                y[i] -= y[j] * R[j,i]
            y[i] /= R[i,i]

        b = array([self.dot(yi, y[-1]) for yi in y[:-1]])
        self.b.append(b)

        for i in range(self.nHomo):
            y[-1] -= y[i] * b[i]

        if _DEBUG_:
            self.y_next.append([yi.copy() for yi in y])

        self.stored_grad.append(pythonCopy(grad))
        self.computed_grad.append(self._compute_grad())

    def grad(self):
        if not _DEBUG_:
            return self._compute_grad()
        else:
            return self.computed_grad[-1]

    def _compute_grad(self):
        window = sin(linspace(0, pi, len(self.R) + 1))**2
        identities = window[:,newaxis,newaxis] * eye(self.nHomo)
        zero = [zeros(self.nHomo)] * (len(self.R) + 1)
        self.a = solveLss(self.R, identities, self.b, zero)

        grad = 0
        win = sin(linspace(0, pi, len(self.a) - 1))**2
        win /= win.mean()
        for i in range(len(self.a) - 1):
            for j in range(self.nHomo):
                grad += win[i] * self.a[i][j] * self.stored_grad[i][j]
            grad += win[i] * self.stored_grad[i][-1]

            if _DEBUG_:
                y_last = self.y_last[i][-1].copy()
                y_next = self.y_next[i][-1].copy()
                for j in range(self.nHomo):
                    y_last += self.y_last[i][j] * self.a[i][j]
                    y_next += self.y_next[i][j] * self.a[i+1][j]
                discontinuity = abs(y_last - y_next).max()
                if discontinuity > 1E-9:
                    print('Discontinuity in y = ', discontinuity)

        return grad

    def lyapunovExponents(self):
        lam = array([log(diag(Ri)) for Ri in self.R])
        return lam.mean(0)

