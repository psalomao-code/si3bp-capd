#!/usr/bin/env python3
from fractions import Fraction as Q

beta=Q(3,100)
alpha=Q(12,97)
w2=Q(25831,1800)
c=1+2*alpha
rL=Q(36,100)
rR=Q(39,100)
pL=Q(85767,10000)
pR=Q(85770,10000)
zmax=Q(1,10000)

def phi_z0(r,p):
    return p*p/2 + w2/(2*r*r) - 1/r - 4/(alpha*r) + 1

# At r=.36, Phi increases with p>0 and |z|, so minimum is p=pL,z=0.
fL=phi_z0(rL,pL)

# sqrt(r^2+c z^2) <= r + c z^2/(2r).  Since -K/sqrt(.) is increasing
# with the denominator, this gives a rational upper bound for Phi at r=.39.
S=rR+c*zmax*zmax/(2*rR)
fR_upper=pR*pR/2+w2/(2*rR*rR)-1/rR-4/(alpha*S)+1

# H_r <= (r/beta-w2)/r^3.  This upper bound is increasing on [.36,.39],
# hence its maximum is attained at r=.39.
d_upper=(rR/beta-w2)/(rR**3)

print('Phi(0.36) uniform lower bound =', fL, '=', float(fL))
print('Phi(0.39) uniform upper bound =', fR_upper, '=', float(fR_upper))
print('H_r uniform upper bound        =', d_upper, '=', float(d_upper))

assert fL > Q(552,1000)
assert fR_upper < -Q(512,1000)
assert d_upper < -Q(227,10)
print('ANALYTIC PRECHECKS PASSED')
