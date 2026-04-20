 PDF To Markdown Converter
Debug View
Result View
1 SUPPLEMENTARY APPENDIX
Differential Equations

Temporal evolution equation for transmembrane voltage in single cell:
dV
dt
=−
Iion+Istim
Cm
(S1)
Total ionic currentIion:
Iion=IN a+IK 1 +IClCa+IKur+IKr+IKs
+ICa,L+Ip,Ca+IN aK
+IN aCa+Ib,N a+Ib,Ca
(S2)
Spatiotemporal evolution equation for transmembrane voltage in tissue:
∂V
∂t
=−
Iion+Istim
Cm
+D∇^2 V (S3)
Temporal evolution equation for all time-dependent gating variables:
dn
dt
=
n∞−n
τn
(S4)
Temporal evolution equations for the different ionic concentrations:
d[Na+]i
dt
=
− 3 IN aK− 3 IN aCa−Ib,N a−IN a
Fvi
·Cm (S5)
d[K+]i
dt
=
2 IN aK−IK 1 −IKur−IKr−IKs−Ib,K
Fvi
·Cm (S6)
d[Cl−]i
dt
=
IClCa
F·vi
·Cm (S7)
d[Ca2+]i
dt
=
B 1
B 2
(S8)
B1 =
2 IN aCa−Ip,Ca−ICa,L−Ib,Ca
2 Fvi
·Cm+
vup(Iup,leak−Iup) +Irelvrel
vi
(S9)
1
B2 = 1 +
[Trpn]maxKm,T rpn
([Ca2+]i+Km,T rpn)^2
+
[Cmdn]maxKm,Cmdn
([Ca2+]i+Km,Cmdn)^2
(S10)
d[Ca2+]up
dt
=Iup−Iup,leak−Itr
vrel
vup
(S11)
d[Ca2+]rel
dt
= (Itr−Irel)
(
1 +
[Csqn]maxKm,Csqn
([Ca2+]rel+Km,Csqn)^2
)− 1
(S12)
Ionic Reversal Potential

EX= (Vi−Vo) =
RT
zF
ln(
[X]o
[X]i
), X=Na+,K+,Ca2+,Cl− (S13)
Ionic Currents

FastNa+Current

IN a=gN am^3 hj(V−EN a) (S14)
αm= 0. 32
V+ 47. 13
1 −exp [− 0 .1(V+ 47.13)]
, if V ̸= 47. 13
αm= 3. 2 , if V =− 47. 13
(S15)
τm=
1
αm+βm
· 1. 7 (S16)
βm= 0.08 exp

−
V
11

(S17)
αh= 0.135 exp

−
V+ 80
6. 8

, if V <− 40
αh= 0, if V ≥− 40
(S18)
βh= 3.56 exp(0. 079 V) + 3. 1 · 105 exp(0. 35 V), if V <− 40
βh=
0. 13
"
1 + exp

−
V+ 10. 66
11. 1

#!− 1
, if V ≥− 40
(S19)
τh=
1
αh+βh
· 2 (S20)
2
αj= [−127140 exp(0. 2444 V)− 3. 474 · 10 −^5 exp(− 0. 04391 V)]
·
V+ 37. 78
1 + exp(0.311(V+ 79.23))
, if V <− 40
αj= 0, if V ≥− 40
(S21)
βj= 0. 1212
exp(− 0. 01052 V)
1 + exp(− 0 .1378(V+ 40.14))
, if V <− 40
βj= 0. 3
exp(− 2. 535 · 10 −^7 )
1 + exp(− 0 .1(V+ 32))
, if V ≥− 40
(S22)
τj=
1
αj+βj
· 2 (S23)
m∞=
1
1 + exp(−(V 7 +43). 7
(S24)
h∞=
1
1 + exp((V+66 4. 4 .5)
(S25)
j∞=
1
1 + exp((V+66 4. 1 .5)
(S26)
Inward Rectifier,IK 1

IK 1 =gK 1 ·
V−EK− 5
1 +exp(0. 063 ·(V+ 70))
(S27)
Transient Outward: Calcium-Driven Chloride Current,Ito=Ito, 2 ≡ICl,Ca

IClCa=gCl,CaqCa(V−ECl) (S28)
qCa,∞= 1−
"
1
1 +

Fn
1. 1 e− 10
 3
(S29)
τCa= 2 (S30)
Ultrarrapid Delayed Rectifier Current,IKur

IKur=gKur·u^3 a·(0. 25 ui,f+ 0. 75 ui,s)·(V−EK) (S31)
Frontiers 3

gKur=gKur,amp
h
0 .005 +
0. 05
1 + exp(−V 13 −^15 )
i
(S32)
αu(a)= 0. 65
"
exp(−
V+ 10
8. 5
) + exp(−
V− 30
59. 0
)
#− 1
(S33)
βu(a)= 0. 65
"
2 .5 + exp(
V+ 82
17. 0
)
#− 1
(S34)
τu(a)=
1
KQ 10 (αu(a)+βu(a))
(S35)
ua(∞)=
1
1 + exp(−V+30 9. 6.^3 )
(S36)
τuif= 400 + 1068e−(
V
50 )
2
(S37)
τuis= 2000 + 60000e−(
V+39. 3
30 )
2
(S38)
ui,f(∞)=ui,f(∞)=
1
1 + exp((V+17 5. 849 .358))
(S39)
Rapid Delayed Rectifier Current,IKr

IKr=gKr·xr·
V−EK
1 +exp(V− 879. 2217.^4825 )
(S40)
αx,r= 0. 0003
V+ 14. 1
1 −exp(−V+14 5.^1 )
(S41)
βx,r= 7. 3898 · 10 −^5
V− 3. 3328
exp(V− 5.^31237.^3328 )− 1
(S42)
τx,s=
1
αx(r)+βx(r)
(S43)
xr,∞=
1
1 + exp(−(V− 9.^433047 .445095))
(S44)
Slow Delayed Rectifier,IKs

IKs=gKs·x^2 s·(V−EK) (S45)
4
αx,s= 4· 10 −^5
V − 18. 80816
1 −exp(−V−^1817.^80816 )
(S46)
βx,s= 3. 5 · 10 −^5
V− 18. 80816
exp(V−^189.^80816 )− 1
(S47)
τx,s=
1
αx(s)+βx(s)
· 0. 5 (S48)
xs,∞=
h 1
1 + exp(−(V− 1218. 6475 .80816))
i 1 / 2
(S49)
L-TypeCa2+Current,ICa,L

ICa,L=gCa,LdffCa(V− 65 .0) (S50)
τd=
1 −exp(−V 6 .+5 24 )
0 .035(V+ 5)
"
1 + exp(−V 6 .+5 24 )
# (S51)
d∞=
1
1 + exp(−(V+5) 8 )
(S52)
τf= 9·[0.0197 exp(− 0. 03372 (V+ 5)^2 ) + 0.02]−^1 (S53)
f∞=
1
1 + exp(V 6 +28. 9 )
(S54)
τf(Ca)= 2, fCa,∞=
1
1 +[Ca
2+]i
0. 00035
(S55)
NaKPump Current

IN aK=IN aK,maxfN aK
1
1 +{Km,N a(i)/[Na+]i}^1.^5
·
[K+]o
[K+]o+Km,K(o)
(S56)
fN aK=
"
1 + 0.1245 exp

− 0. 1
FV
RT

+ 0. 0365 σexp

−
FV
RT

#− 1
(S57)
σ=
1
7
"
exp
[Na+]
o
67. 3

− 1
(S58)
Frontiers 5

Na+/Ca2+Exchanger

IN aCa=IN aCa,max·
B 1
B 2
(S59)
B1 = exp[γV F/(RT)][Na+]^3 i[Ca2+]o−exp[(γ−1)V F/(RT)][Na+]^3 o[Ca2+]i (S60)
B2 = (Km,N a^3 + [Na+]^3 o)(Km,Ca+ [Ca2+]o)·{1 +ksatexp[(γ−1)V F/(RT)]} (S61)
Background Currents

Ib,X=gb,X(V−EX), X=Na+,Ca2+ (S62)
Ca2+Pump Current

Ip,Ca=Ip,Ca(max)
[Ca2+]i
0 .0005 + [Ca2+]i
(S63)
Ca2+Release Current from the JSR,Irel

Irel=krel·u^2 vw·([Ca2+]rel−[Ca2+]i) (S64)
τu= 8, u∞=
1
1 + exp

Fn− 1. 367 · 10 −^13
13. 67 · 10 −^16
 (S65)
τv= 1.91 + 2. 09
1 + exp
"
−
Fn− 1. 367 · 10 −^13
13. 67 · 10 −^16
#!− 1
(S66)
v∞= 1−
1
1 + exp

−Fn−^6.^835 ·^10
− 14
13. 67 · 10 −^16
 (S67)
τw= 6. 0 ·
1 −exp(−V− 57.^9 )
"
1 + 0.3 exp(−V− 57.^9 )
(V− 7 .9)
(S68)
w∞=
1
1 + exp(V− 870 )
(S69)
Fn= 10−^12 vrelIrel−
5 · 10 −^13
F
(
1
2
ICa,L−
1
5
IN aCa)·Cm (S70)
6
Transfer Current from NSR to JSR

Itr=
[Ca2+]up−[Ca2+]rel
τtr
(S71)
τtr= 180 (S72)
Ca2+Uptake by the NSR

Iup=
Iup(max)
1 + (Kup/[Ca2+]i)
(S73)
Ca2+Leak Current by the NSR

Iup,leak=
[Ca2+]up
[Ca2+]up(max)
Iup(max) (S74)
Calcium Buffers

[Ca2+]Cmdn= [Cmdn]max
[Ca2+]i
[Ca2+]i+Km,Cmdn
(S75)
[Ca2+]T rpn= [Trpn]max
[Ca2+]i
[Ca2+]i+Km,T rpn
(S76)
[Ca2+]Csqn= [Csqn]max
[Ca2+]rel
[Ca2+]rel+Km,Csqn
(S77)
Numerical Integration

For the gating variables, at timet:
n(t+1)=n∞−[n∞−n(t)]e−
∆t
τ (S78)
For the rest of differential equations:

f(t+1)=f(t)+ ∆t
∂f
∂t
(S79)
Frontiers 7

2 SUPPLEMENTARY TABLES
Table S1.Set of model parameters and universal constants

Parameter Defintion Value
R Gas constant 8.3143JK−^1 mol−^1
T Temperature 310 K
F Faraday constant 96. 4867 C/mmol
Cm Membrane capacitance 100 pF
D Diffusion Coefficient 0. 00126 cm^2 /ms
vcell Cell volume 20100 μm^3
vi Intracellular volume 13668 μm^3
vup SR uptake compartment volume 1109. 52 μm^3
vrel SR release compartment volume 96. 48 μm^3
[K+]o ExtracellularK+concentration 5. 4 mM
[Na+]o ExtracellularN a+concentration 140 mM
[Ca2+]o ExtracellularCa2+concentration 1. 8 mM
[Cl−]o ExtracellularCl−concentration 132 mM
gN a MaximalINaconductance 13. 9900 nS/pF
gK 1 MaximalIK 1 conductance 0. 08218 nS/pF
gto, 2 /gCl,Ca MaximalItoconductance 0. 15731 nS/pF
gKur,amp MaximalIKurconductance 0. 45539 nS/pF
gKr MaximalIKrconductance 0. 01730 nS/pF
gKs MaximalIKsconductance 0. 0594 nS/pF
gCa,L MaximalICa,Lconductance 0. 06574 nS/pF
gb,Ca MaximalIb,Caconductance 0. 00113 nS/pF
gb,N a MaximalIb,Naconductance 0. 000674 nS/pF
IN aK(max) MaximalINaK 0. 94935 pA/pF
IN aCa(max) MaximalINaCa 2304 pA/pF
Ip,Ca(max) MaximalIp,Ca 0. 275 pA/pF
Iup(max) MaximalIup 0. 005 mM/ms
KQ 10 Temperature scaling factor 3
γ Voltage dependence parameter forINaCa 0.
Km,N a(i) [N a+]ihalf-saturation constant forINaK 10 mM
Km,K(o) [K+]ohalf-saturation constant forINaK 1. 5 mM
Km,N a [N a+]ohalf-saturation constant forINaCa 87. 5 mM
Km,Ca [Ca2+]ohalf-saturation constant forINaCa 1. 38 mM
ksat Saturation factor forINaCa 0.
krel Maximal release rate forIrel 30 ms−^1
Kup [Ca2+]ihalf-saturation constant forIup 0. 00092 mM
[Ca2+]up(max) MaximalCa2+concentration in uptake compartment 15 mM
[Cmdn]max Total calmodulin concentration in myoplasm 0. 05 mM
[Trpn]max Total troponin concentration in myoplasm 0. 07 mM
[Csqn]max Total calsequestrin concentration in SR release compartment 10 mM
Km,Cmdn [Ca2+]ihalf-saturation constant for calmodulin 0. 00238 mM
Km,T rpn [Ca2+]ihalf-saturation constant for troponin 0. 0005 mM
Km,Csqn [Ca2+]ihalf-saturation constant forIup 0. 8 mM
8
Table S2.Set of initial values of variables used in the model.

Variable Value
V -76.
h 0 9. 65 × 10 −^1
d 0 1. 37 × 10 −^4
xr, 0 3. 29 × 10 −^5
[Na+]i 1. 117 × 10 −^1
[K+]i 1. 39 × 102
[Ca2+]rel 1. 488
oi, 0 9. 99 × 10 −^1
uif, 0 1
uis, 0 1
[Cmdn]i 2. 05 × 10 −^3
[Csqn]i 6. 51
[Trpn]i 1. 18 × 10 −^2
v 0 1
m 0 2. 91 × 10 −^3
j 0 9. 78 × 10 −^1
f 0 9. 99 × 10 −^1
xs, 0 1. 87 × 10 −^2
[Ca2+]i 0. 0001013
[Ca2+]up 1. 488
[Cl−]i 30
oa, 0 3. 04 × 10 −^2
ua, 0 4. 96 × 10 −^3
fCa, 0 7. 75 × 10 −^1
u 0 2. 35 × 10 −^2
w 0 9. 99 × 10 −^1
Irel 0
Ib,k 0
qCa, 0 0.
This is a offline tool, your data stays locally and is not send to any server!
Feedback & Bug Reports