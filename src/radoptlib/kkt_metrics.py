"""Check a candidate solution of  minimize 0.5 x'Qx + q'x + c  subject to  Gx <= h.

All residuals are relative and use the infinity norm:
  primal_res = ||max(Gx - h, 0)|| / (1 + max(||Gx||, ||h||))
  dual_res   = ||Qx + q + G'y|| / (1 + max(||Qx||, ||q||, ||G'y||))
  gap        = |pobj - dobj| / (1 + max(|pobj|, |dobj|)),
               pobj = 0.5 x'Qx + q'x + c,  dobj = -0.5 x'Qx - h'y + c
where y >= 0 are the multipliers of Gx <= h (Lagrangian 0.5 x'Qx + q'x + c + y'(Gx - h)).
"""
import numpy as np


def _inf(v):
    return float(np.linalg.norm(v, np.inf)) if v.size else 0.0


def verify_solution(problem, x, dual=None, tol=None):
    """Residuals of a primal point x and (optionally) dual y >= 0 for Gx <= h.

    Returns a dict with objective, primal_res and, when dual is given, dual_res, gap and dual_infeas
    (= max(-y, 0), how far y is from being nonnegative). With tol, adds 'ok' = every reported
    residual <= tol."""
    x = np.asarray(x, dtype=np.float64).ravel()
    if x.shape[0] != problem.n:
        raise ValueError('x has length {}, problem has n = {}'.format(x.shape[0], problem.n))
    Q, q, c, G, h = problem.Q, problem.q, problem.c, problem.G, problem.h
    Qx, Gx = Q @ x, G @ x
    pobj = float(0.5 * x @ Qx + q @ x + c)
    out = {'objective': pobj,
           'primal_res': _inf(np.maximum(Gx - h, 0.0)) / (1.0 + max(_inf(Gx), _inf(h)))}
    if dual is not None:
        y = np.asarray(dual, dtype=np.float64).ravel()
        if y.shape[0] != problem.m:
            raise ValueError('dual has length {}, problem has m = {}'.format(y.shape[0], problem.m))
        Gty = G.T @ y
        dobj = float(-0.5 * x @ Qx - h @ y + c)
        out['dual_res'] = _inf(Qx + q + Gty) / (1.0 + max(_inf(Qx), _inf(q), _inf(Gty)))
        out['gap'] = abs(pobj - dobj) / (1.0 + max(abs(pobj), abs(dobj)))
        out['dual_infeas'] = _inf(np.maximum(-y, 0.0))
    if tol is not None:
        out['ok'] = all(out[k] <= tol for k in ('primal_res', 'dual_res', 'gap', 'dual_infeas') if k in out)
    return out