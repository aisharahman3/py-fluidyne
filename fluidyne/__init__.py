import math

def reynolds(rho, v, d, mu):
    return rho * v * d / mu

def is_laminar(re):
    return re < 2300.0

def friction_factor(re, eps, d, tol=1e-8):
    """Darcy friction factor.  Laminar: 64/Re; turbulent: Colebrook-White
    solved by fixed-point iteration."""
    if re < 2300:
        return 64.0 / re
    rr = eps / d
    f = 0.02
    for _ in range(100):
        rhs = -2.0 * math.log10(rr / 3.7 + 2.51 / (re * math.sqrt(f)))
        fn = 1.0 / rhs ** 2
        if abs(fn - f) < tol:
            return fn
        f = fn
    return f

def head_loss(f, length, d, v, g=9.80665):
    return f * (length / d) * v * v / (2 * g)

def bernoulli(p1, rho, v1, z1, v2, z2, g=9.80665):
    return p1 + 0.5 * rho * (v1 * v1 - v2 * v2) + rho * g * (z1 - z2)

def hagen_poiseuille(dp, radius, length, mu):
    return math.pi * radius ** 4 * dp / (8.0 * mu * length)
