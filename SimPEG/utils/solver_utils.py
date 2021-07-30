from __future__ import print_function
import numpy as np
from scipy.sparse import linalg
from .mat_utils import mkvc
import warnings
import inspect


def _checkAccuracy(A, b, X, accuracyTol):
    nrm = np.linalg.norm(mkvc(A * X - b), np.inf)
    nrm_b = np.linalg.norm(mkvc(b), np.inf)
    if nrm_b > 0:
        nrm /= nrm_b
    if nrm > accuracyTol:
        msg = "### SolverWarning ###: Accuracy on solve is above tolerance: {0:e} > {1:e}".format(
            nrm, accuracyTol
        )
        print(msg)
        warnings.warn(msg, RuntimeWarning)


def SolverWrapD(fun, factorize=True, checkAccuracy=True, accuracyTol=1e-6, name=None):
    """
    Wraps a direct Solver.

    ::

        import scipy.sparse as sp
        Solver   = solver_utils.SolverWrapD(sp.linalg.spsolve, factorize=False)
        SolverLU = solver_utils.SolverWrapD(sp.linalg.splu, factorize=True)

    """

    def __init__(self, A, **kwargs):
        self.A = A.tocsc()

        self.checkAccuracy = kwargs.pop("checkAccuracy", checkAccuracy)
        self.accuracyTol = kwargs.pop("accuracyTol", accuracyTol)

        func_params = inspect.signature(fun).parameters
        # First test if function excepts **kwargs,
        # in which case we do not need to cull the kwargs
        do_cull = True
        for param_name in func_params:
            param = func_params[param_name]
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                do_cull = False
        if do_cull:
            # build a dictionary of valid kwargs
            culled_args = {}
            for item in kwargs:
                if item in func_params:
                    culled_args[item] = kwargs[item]
                else:
                    warnings.warn(
                        f"{item} is not a valid keyword for {fun.__name__} and will be ignored"
                    )
            kwargs = culled_args

        self.kwargs = kwargs

        if factorize:
            self.solver = fun(self.A, **kwargs)

    def __mul__(self, b):
        if not isinstance(b, np.ndarray):
            raise TypeError("Can only multiply by a numpy array.")

        if len(b.shape) == 1 or b.shape[1] == 1:
            b = b.flatten()
            # Just one RHS

            if b.dtype is np.dtype("O"):
                b = b.astype(type(b[0]))

            if factorize:
                X = self.solver.solve(b, **self.kwargs)
            else:
                X = fun(self.A, b, **self.kwargs)
        else:  # Multiple RHSs
            if b.dtype is np.dtype("O"):
                b = b.astype(type(b[0, 0]))

            X = np.empty_like(b)

            for i in range(b.shape[1]):
                if factorize:
                    X[:, i] = self.solver.solve(b[:, i])
                else:
                    X[:, i] = fun(self.A, b[:, i], **self.kwargs)

        if self.checkAccuracy:
            _checkAccuracy(self.A, b, X, self.accuracyTol)
        return X

    def __matmul__(self, other):
        return self * other

    def clean(self):
        if factorize and hasattr(self.solver, "clean"):
            return self.solver.clean()

    return type(
        name if name is not None else fun.__name__,
        (object,),
        {
            "__init__": __init__,
            "clean": clean,
            "__mul__": __mul__,
            "__matmul__": __matmul__,
        },
    )


def SolverWrapI(fun, checkAccuracy=True, accuracyTol=1e-5, name=None):
    """
    Wraps an iterative Solver.

    ::

        import scipy.sparse as sp
        SolverCG = solver_utils.SolverWrapI(sp.linalg.cg)

    """

    def __init__(self, A, **kwargs):
        self.A = A

        self.checkAccuracy = kwargs.pop("checkAccuracy", checkAccuracy)
        self.accuracyTol = kwargs.pop("accuracyTol", accuracyTol)

        func_params = inspect.signature(fun).parameters
        # First test if function excepts **kwargs,
        # in which case we do not need to cull the kwargs
        do_cull = True
        for param_name in func_params:
            param = func_params[param_name]
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                do_cull = False
        if do_cull:
            # build a dictionary of valid kwargs
            culled_args = {}
            for item in kwargs:
                if item in func_params:
                    culled_args[item] = kwargs[item]
                else:
                    warnings.warn(
                        f"{item} is not a valid keyword for {fun.__name__} and will be ignored"
                    )
            kwargs = culled_args

        self.kwargs = kwargs

    def __mul__(self, b):
        if not isinstance(b, np.ndarray):
            raise TypeError("Can only multiply by a numpy array.")

        if len(b.shape) == 1 or b.shape[1] == 1:
            b = b.flatten()
            # Just one RHS
            out = fun(self.A, b, **self.kwargs)
            if isinstance(out, tuple) and len(out) == 2:
                # We are dealing with scipy output with an info!
                X = out[0]
                self.info = out[1]
            else:
                X = out
        else:  # Multiple RHSs
            X = np.empty_like(b)
            for i in range(b.shape[1]):
                out = fun(self.A, b[:, i], **self.kwargs)
                if isinstance(out, tuple) and len(out) == 2:
                    # We are dealing with scipy output with an info!
                    X[:, i] = out[0]
                    self.info = out[1]
                else:
                    X[:, i] = out

        if self.checkAccuracy:
            _checkAccuracy(self.A, b, X, self.accuracyTol)
        return X

    def __matmul__(self, other):
        return self * other

    def clean(self):
        pass

    return type(
        name if name is not None else fun.__name__,
        (object,),
        {
            "__init__": __init__,
            "clean": clean,
            "__mul__": __mul__,
            "__matmul__": __matmul__,
        },
    )


Solver = SolverWrapD(linalg.spsolve, factorize=False, name="Solver")
SolverLU = SolverWrapD(linalg.splu, factorize=True, name="SolverLU")
SolverCG = SolverWrapI(linalg.cg, name="SolverCG")
SolverBiCG = SolverWrapI(linalg.bicgstab, name="SolverBiCG")


class SolverDiag(object):
    """docstring for SolverDiag"""

    def __init__(self, A, **kwargs):
        self.A = A
        self._diagonal = A.diagonal()
        for kwarg in kwargs:
            warnings.warn(f"{kwarg} is not recognized and will be ignored")

    def __mul__(self, rhs):
        n = self.A.shape[0]
        assert rhs.size % n == 0, "Incorrect shape of rhs."
        nrhs = rhs.size // n

        if len(rhs.shape) == 1 or rhs.shape[1] == 1:
            x = self._solve1(rhs)
        else:
            x = self._solveM(rhs)

        if nrhs == 1:
            return x.flatten()
        elif nrhs > 1:
            return x.reshape((n, nrhs), order="F")

    def __matmul__(self, other):
        return self * other

    def _solve1(self, rhs):
        return rhs.flatten() / self._diagonal

    def _solveM(self, rhs):
        n = self.A.shape[0]
        nrhs = rhs.size // n
        return rhs / self._diagonal.repeat(nrhs).reshape((n, nrhs))

    def clean(self):
        pass
