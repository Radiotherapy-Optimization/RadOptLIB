"""Check a candidate solution of  minimize 0.5 x'Qx + q'x + c  subject to  Gx <= h,  l <= x <= u.

Multipliers (all >= 0): dual_ineq_cons (y) for Gx <= h, dual_lower_bound (z_l) for x >= l,
dual_upper_bound (z_u) for x <= u.
Lagrangian 0.5 x'Qx + q'x + c + y'(Gx - h) + z_l'(l - x) + z_u'(x - u). A bound that is missing
(l_j = -inf or u_j = +inf) must have a zero multiplier.

All residuals are relative and use the infinity norm:
  primal_res = max(||max(Gx - h, 0)||, ||max(l - x, 0)||, ||max(x - u, 0)||)
               / (1 + max(||Gx||, ||h||, ||x||, ||l||, ||u||))          (finite bounds only)
  dual_res   = ||Qx + q + G'y - z_l + z_u|| / (1 + max(||Qx||, ||q||, ||G'y||, ||z_l||, ||z_u||))
  gap        = |pobj - dobj| / (1 + max(|pobj|, |dobj|)),
               pobj = 0.5 x'Qx + q'x + c,  dobj = -0.5 x'Qx - h'y + l'z_l - u'z_u + c
  dual_infeas = how far the multipliers are from being nonnegative (and zero at missing bounds).
"""
import numpy as np


def _inf(v):
    return float(np.linalg.norm(v, np.inf)) if v.size else 0.0


def _vec(v, size, what):
    v = np.zeros(size) if v is None else np.asarray(v, dtype=np.float64).ravel()
    if v.shape[0] != size:
        raise ValueError('{} has length {}, expected {}'.format(what, v.shape[0], size))
    return v


def verify_solution(problem, x, dual_ineq_cons=None, dual_lower_bound=None, dual_upper_bound=None, tol=None):
    """Residuals of a primal point x and (optionally) multipliers.

    Returns a dict with objective and primal_res and, when any multiplier is given, dual_res, gap and
    dual_infeas (multipliers not given count as zero). With tol, adds 'ok' = every reported
    residual <= tol. Example with a reference solution:
        sol = load_solution(name)
        verify_solution(load(name), sol['x'], sol['dual_ineq_cons'], sol['dual_lower_bound'],
                        sol['dual_upper_bound'])"""
    n, m = problem.n, problem.m
    x = _vec(x, n, 'x')
    Q, q, c, G, h = problem.Q, problem.q, problem.c, problem.G, problem.h
    l = np.asarray(getattr(problem, 'l', np.full(n, -np.inf)), dtype=np.float64)
    u = np.asarray(getattr(problem, 'u', np.full(n, np.inf)), dtype=np.float64)
    fl, fu = np.isfinite(l), np.isfinite(u)
    Qx, Gx = Q @ x, G @ x
    pobj = float(0.5 * x @ Qx + q @ x + c)
    viol = max(_inf(np.maximum(Gx - h, 0.0)), _inf(np.maximum(l[fl] - x[fl], 0.0)), _inf(np.maximum(x[fu] - u[fu], 0.0)))
    out = {'objective': pobj,
           'primal_res': viol / (1.0 + max(_inf(Gx), _inf(h), _inf(x), _inf(l[fl]), _inf(u[fu])))}
    if not (dual_ineq_cons is None and dual_lower_bound is None and dual_upper_bound is None):
        y = _vec(dual_ineq_cons, m, 'dual_ineq_cons')
        zl = _vec(dual_lower_bound, n, 'dual_lower_bound')
        zu = _vec(dual_upper_bound, n, 'dual_upper_bound')
        Gty = G.T @ y
        dobj = float(-0.5 * x @ Qx - h @ y + l[fl] @ zl[fl] - u[fu] @ zu[fu] + c)
        out['dual_res'] = _inf(Qx + q + Gty - zl + zu) / (1.0 + max(_inf(Qx), _inf(q), _inf(Gty), _inf(zl), _inf(zu)))
        out['gap'] = abs(pobj - dobj) / (1.0 + max(abs(pobj), abs(dobj)))
        out['dual_infeas'] = max(_inf(np.maximum(-y, 0.0)), _inf(np.maximum(-zl, 0.0)), _inf(np.maximum(-zu, 0.0)),
                                 _inf(zl[~fl]), _inf(zu[~fu]))
    if tol is not None:
        out['ok'] = all(out[k] <= tol for k in ('primal_res', 'dual_res', 'gap', 'dual_infeas') if k in out)
    return out