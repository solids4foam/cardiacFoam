import sympy as sp
import re, sys, math, random

X, Y, Z, t = sp.symbols('X Y Z t', real=True)
Ax, Ay, Az = sp.symbols('Ax Ay Az', real=True)
Tmax, V0, gamma = sp.symbols('Tmax V0 gamma', real=True)
mu, K, rho = sp.symbols('mu K rho', real=True)
pi = sp.pi
Xv = [X, Y, Z]

st = sp.sin(t)
D = sp.Matrix([Ax*X**2*Y*st, Ay*Y**2*Z*st, Az*Z**2*X*st])

# F = I + Grad_X D   (F_ij = dD_i/dX_j)
F = sp.eye(3) + sp.Matrix(3, 3, lambda i, j: sp.diff(D[i], Xv[j]))
J = F.det()

# --- passive compressible neo-Hookean, EXACTLY as solids4foam neoHookeanElastic ---
b = F*F.T                                  # left Cauchy-Green
bbar = J**sp.Rational(-2, 3) * b
dev = lambda A: A - sp.trace(A)/3*sp.eye(3)
s = mu*dev(bbar)
sigma_hyd = sp.Rational(1, 2)*K*(J**2 - 1)
sigma_p = (sigma_hyd*sp.eye(3) + s)/J

# --- active fibre stress, f0 = (1,0,0) ---
f0 = sp.Matrix([1, 0, 0])
Vm = sp.sqrt(1 + t)*sp.cos(pi*X)*sp.cos(2*pi*Y)*sp.cos(3*pi*Z)
lam = sp.sqrt((F*f0).dot(F*f0))
Ta = Tmax*Vm**2/(V0**2 + Vm**2)*(1 + gamma*(lam - 1))
sigma_a = (F*(Ta*(f0*f0.T))*F.T)/J

sigma = sigma_p + sigma_a
P = J*sigma*F.inv().T                       # 1st Piola-Kirchhoff

DivP = sp.Matrix([sum(sp.diff(P[i, j], Xv[j]) for j in range(3)) for i in range(3)])
inertia = rho*sp.diff(D, t, 2)              # = -rho*D
B = inertia - DivP                          # B = rho*d2D/dt2 - Div(P)   (g=0)

args = (X, Y, Z, t, Ax, Ay, Az, Tmax, V0, gamma, mu, K, rho)
fB = sp.lambdify(args, [B[0], B[1], B[2]], 'math')
print("sympy B built OK")

# Emit sympy's C++ form (for optional install) using cxxcode/ccode
from sympy.printing.cxx import cxxcode
with open('/tmp/B_expr_sympy.H', 'w') as fh:
    for nm, expr in (('Bx', B[0]), ('By', B[1]), ('Bz', B[2])):
        fh.write("const scalar %s = %s;\n\n" % (nm, cxxcode(expr, standard='c++11')))
print("wrote /tmp/B_expr_sympy.H")

# Fixed test points; print sympy B for the C++ harness to match
PTS = [(0.2, 0.3, 0.4, 0.1), (0.55, 0.65, 0.45, 0.05), (0.8, 0.15, 0.7, 0.08)]
base = dict(Ax=0.02, Ay=0.02, Az=0.02, V0=1.0, gamma=1.0, mu=3846.15, K=8333.33, rho=1060.0)
for label, Tm in (("PASSIVE_Tmax0", 0.0), ("FULL_Tmax1000", 1000.0)):
    for (xx, yy, zz, tt) in PTS:
        v = fB(xx, yy, zz, tt, base['Ax'], base['Ay'], base['Az'], Tm,
               base['V0'], base['gamma'], base['mu'], base['K'], base['rho'])
        print("SYMPY %s X=%.2f Y=%.2f Z=%.2f t=%.2f  B=(% .10e % .10e % .10e)"
              % (label, xx, yy, zz, tt, v[0], v[1], v[2]))
import sys; sys.exit(0)

# ---- load committed B_expr.H as an evaluable python function ----
def load_committed(path):
    txt = open(path).read().replace('std::', '')
    exprs = {}
    for name in ('Bx', 'By', 'Bz'):
        m = re.search(r'const scalar %s\s*=\s*(.*?);' % name, txt, re.S)
        exprs[name] = m.group(1)
    ns = {'sin': math.sin, 'cos': math.cos, 'sqrt': math.sqrt,
          'pow': math.pow, 'M_PI': math.pi}
    def f(X_, Y_, Z_, t_, Ax_, Ay_, Az_, Tmax_, V0_, gamma_, mu_, K_, rho_):
        loc = dict(X=X_, Y=Y_, Z=Z_, t=t_, Ax=Ax_, Ay=Ay_, Az=Az_,
                   Tmax=Tmax_, V0=V0_, gamma=gamma_, mu=mu_, K=K_, rho=rho_)
        loc.update(ns)
        return [eval(exprs[n], {'__builtins__': {}}, loc) for n in ('Bx', 'By', 'Bz')]
    return f

committed = sys.argv[1] if len(sys.argv) > 1 else None
fC = load_committed(committed) if committed else None

def relerr(a, b):
    d = abs(a-b); s = max(abs(a), abs(b), 1e-30)
    return d/s

random.seed(1)
for label, Tm in (("PASSIVE (Tmax=0)", 0.0), ("FULL (Tmax=1000)", 1000.0)):
    print("\n===", label, "===")
    worst = 0.0
    for _ in range(6):
        vals = dict(
            X_=random.uniform(0.05, 0.95), Y_=random.uniform(0.05, 0.95),
            Z_=random.uniform(0.05, 0.95), t_=random.uniform(0.01, 0.1),
            Ax_=0.02, Ay_=0.02, Az_=0.02, Tmax_=Tm, V0_=1.0, gamma_=1.0,
            mu_=3846.15, K_=8333.33, rho_=1060.0)
        a = fB(vals['X_'], vals['Y_'], vals['Z_'], vals['t_'], vals['Ax_'], vals['Ay_'],
               vals['Az_'], vals['Tmax_'], vals['V0_'], vals['gamma_'], vals['mu_'], vals['K_'], vals['rho_'])
        line = "  pt B_sympy=(%+.4g %+.4g %+.4g)" % tuple(a)
        if fC:
            c = fC(**vals)
            re3 = [relerr(a[i], c[i]) for i in range(3)]
            worst = max(worst, max(re3))
            line += "  committed=(%+.4g %+.4g %+.4g)  relerr=(%.2e %.2e %.2e)" % (c[0], c[1], c[2], re3[0], re3[1], re3[2])
        print(line)
    if fC:
        print("  -> worst rel err this block: %.3e" % worst)
