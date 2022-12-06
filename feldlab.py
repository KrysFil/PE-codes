import sympy
import numpy as np
from matplotlib import pyplot as plt
from sympy import symbols, simplify, lambdify, I, pi, Abs, expand, log

f,R,C,L,Rc,Cc = symbols("f R C L Rc Cc", real=True)
f_r,Q,f_0 = symbols("f_r Q f_0", real=True)

freq = np.logspace(1,4,1000)

RL_C = 1/(1+I*2*pi*R*C-(2*pi*f)**2*L*C)
C_Rc = (I * 2 * pi * f * Rc * Cc) / (1 + I * 2 * pi * f * Rc * Cc)


funk = RL_C*C_Rc
sub = funk.subs([(R*C,1/(2*pi*Q*f_r)),(L*C,1/(2*pi*f_r)**2),
           (Rc*Cc,1/(2*pi*f_0))])
print(sub)
sub = 20*log(Abs(sub),10)
print(sub)
fuss = lambdify([f,f_r,Q,f_0],sub)
plt.plot(freq,fuss(freq,100,10,1000))
plt.semilogx()
plt.show()