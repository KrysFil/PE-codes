import numpy as np
from sympy import symbols, simplify, lambdify, I, pi, Abs, log

# Symbols in sympy functions
f, R, C, L = symbols("f R C L", real=True)
f_0, f_r, Q = symbols("f_0 f_r Q", real=True)

fc, Rc, Cc, Lc = symbols("fc Rc Cc Lc", real=True)
f_0c, f_rc, Qc = symbols("f_0c f_rc Qc", real=True)

# Collection of basic transfer functions
# Classification: __ means parallel, _ is the U_out. In the order of components.

# Low pass
R_C = 1 / (1 + I * 2 * pi * f * R * C)
R_Cc = 1 / (1 + I * 2 * pi * f * Rc * Cc)

# High pass
C_R = (I * 2 * pi * f * R * C) / (1 + I * 2 * pi * f * R * C)
C_Rc = (I * 2 * pi * f * Rc * Cc) / (1 + I * 2 * pi * f * Rc * Cc)

# Low pass resonant
RL_C = 1 / (1 + I * 2 * pi * f * R * C - (2 * pi * f) ** 2 * L * C)
RL_Cc = 1 / (1 + I * 2 * pi * f * Rc * Cc - (2 * pi * f) ** 2 * Lc * Cc)

# High pass resonant
RC_L = (-1 * (2 * pi * f)**2 * L * C) / (1 + I * (2 * pi * f) * R * C - (2 * pi * f)**2 * L * C)
RC_Lc = (-1 * (2 * pi * f)**2 * Lc * Cc) / (1 + I * (2 * pi * f) * Rc * Cc - (2 * pi * f)**2 * Lc * Cc)

# Band stop
R_CL = (1 - (2 * pi * f) ** 2 * L * C) / (1 + I * (2 * pi * f) * R * C - ((2 * pi * f) ** 2) * L * C)
R_CLc = (1 - (2 * pi * f) ** 2 * Lc * Cc) / (1 + I * (2 * pi * f) * Rc * Cc - ((2 * pi * f) ** 2) * Lc * Cc)

# Band pass
LC_R = (I * 2 * pi * f * R * C) / (1 + I * 2 * pi * f * R * C - (2 * pi * f) ** 2 * L * C)
LC_Rc = (I * 2 * pi * f * Rc * Cc) / (1 + I * 2 * pi * f * Rc * Cc - (2 * pi * f) ** 2 * Lc * Cc)

def model(function, **kwargs):
    """Function to convert a sympy function to lambda function for SpecClass. Fill here R,C and L."""
    if len(kwargs) == 2:
        inter = function.evalf(subs={R: kwargs["R"], C: kwargs["C"]})
    elif len(kwargs) == 3:
        inter = function.evalf(subs={R: kwargs["R"], C: kwargs["C"], L: kwargs["L"]})
    elif len(kwargs) == 4:
        inter = function.evalf(subs={R: kwargs["R"], C: kwargs["C"], Rc: kwargs["Rc"], Cc: kwargs["Cc"]})
    elif len(kwargs) == 5:
        inter = function.evalf(subs={R: kwargs["R"], C: kwargs["C"], L: kwargs["L"], Rc: kwargs["Rc"], Cc: kwargs["Cc"]})
    elif len(kwargs) == 6:
        inter = function.evalf(subs={R: kwargs["R"], C: kwargs["C"], L: kwargs["L"],Rc: kwargs["Rc"], Cc: kwargs["Cc"], Lc: kwargs["Lc"]})
    else:
        raise TypeError("Set R,C and ev L, Rc, Cc and Lc")
    return lambdify(f, inter)


def fitModel(function, bode=True):
    """
    Note that second order filter must be normal and copy a first order.
    Substitutions: RC = 1 / 2*pi*f_0
    LC = 1 / (2*pi*f_r)**2
    RC = 1 / Q*f_r
    """
    if len(function.free_symbols) == 3:
        inter = function.subs({R * C: (1 / (2 * pi * f_0))})
        if bode:
            absolute = 20 * log(Abs(inter), 10)
        else:
            absolute = Abs(inter)
        simp = simplify(absolute)
        return lambdify((f, f_0), simp)

    elif len(function.free_symbols) == 5:
        inter = function.subs([(R * C, 1 / (2 * pi * f_0)), (Rc * Cc, 1 / (2 * pi * f_0c))])
        if bode:
            absolute = 20 * log(Abs(inter), 10)
        else:
            absolute = Abs(inter)
        simp = simplify(absolute)
        print(simp)
        return lambdify((f, f_0, f_0c), simp)

    elif len(function.free_symbols) == 4:
        inter = function.subs({R * C: (1 / (Q*f_r*2*pi)), L * C: (1 / (2*np.pi*f_r) ** 2)})
        if bode:
            absolute = 20 * log(Abs(inter), 10)
        else:
            absolute = Abs(inter)
        simp = simplify(absolute)
        return lambdify((f, f_r, Q), simp)

    elif len(function.free_symbols) == 6:
        """Mixing orders"""
        sub = function.subs([(R * C, 1 / (2 * pi * Q * f_r)), (L * C, 1 / (2 * pi * f_r) ** 2),
                             (Rc * Cc, 1 / (2 * pi * f_0))])
        if bode:
            sub = 20 * log(Abs(sub), 10)
        else:
            sub = Abs(sub)
        sub = simplify(sub)
        return lambdify([f, f_r, Q, f_0], sub)

    elif len(function.free_symbols) == 7:
        function = simplify(function)
        inter = function.subs({R * C: (1 / (Q * f_r * 2 * pi)), L * C: (1 / (2 * np.pi * f_r) ** 2),
                               Rc * Cc: (1 / (Qc * f_rc * 2 * pi)), Lc * Cc: (1 / (2 * np.pi * f_rc) ** 2)})
        if bode:
            absolute = 20 * log(Abs(inter), 10)
        else:
            absolute = Abs(inter)
        simp = simplify(absolute)
        return lambdify((f, f_r, Q, f_rc, Qc), simp)


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    freqs = np.logspace(4, 5, 100)

    r = 1e2
    c = 1e-7
    l = 1e-2

    rc = 1e1
    cc = 1e-9
    lc = 1e-3

    f_00 = 1 / (2 * np.pi * r * c)
    f_00c = 1 / (2 * np.pi * rc * cc)

    f_rr = 1 / (2 * np.pi * (l * c) ** 0.5)
    Qq = 2 * np.pi * f_rr * l / r

    f_rrc = 1 / (2 * np.pi * (lc * cc) ** 0.5)
    Qqc = 2 * np.pi * f_rrc * lc / rc

    multi = RL_C*C_Rc
    spekModel = model(multi, R=r, C=c, L=l, Rc=rc, Cc=cc)
    spekFit = fitModel(multi)(freqs, f_rr, Qq, f_00c)

    #spekModel = model(multi, R=r, C=c, L=l, Rc=rc, Cc=cc, Lc=lc)
    #spekFit = fitModel(multi)(freqs, f_rr, Qq, f_rrc, Qqc)

    plt.plot(freqs, 20 * np.log10(np.abs(spekModel(freqs))), c="r")
    plt.plot(freqs, spekFit, c="b")
    plt.semilogx()
    plt.show()