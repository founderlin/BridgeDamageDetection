# input variables
# --------------
nm = 2; # quantity of monitoring sections
rn = 10; # sample quantity in each intervalS
n = 6; # quantity of elements
N = rn^n;
# ParaM = [ 0.21 0.42 0.19 0.38 0.64 0.92 ]
# ParaM = [ 0.2134 0.4257 0.1944 0.3875 0.6429 0.9273 ]
# ParaM = [ 0.2134671 0.4257344 0.1944257 0.3875839 0.6429525 0.9273182 ]
ParaM = [ 0.9 0.1944257 0.4257344 0.9 0.6429525 0.9]

# material parameters
A = 6.5934; # [m2]
I = 4.639; # Iy [m4]
E = 34077*10^6; # [N/m^2]
L = 16; # beam length [m]
F = 10^6; # load [N]
q = 164835; # [N/m]
ea = E*A;
ei = E*I;
el = L/n; # element length
ma_para = [ ea ei el]

# load variant 1
f1 = zeros(3*(n+1), 1);
f1[5,1] = F;

# load variant 2
f2 = zeros(3*(n+1), 1);
f2[8,1] = F;

# load variant 3
f_q = fem.loadVct(q, el, n)
f3 = f_q + f1

# load variant 3
f4 = f_q + f2

f5 = zeros(3*(n+1), 1);
f5[14,1] = F;

f6 = zeros(3*(n+1), 1);
f6[17,1] = F;

f7 = f_q + f5
f8 = f_q + f6
f9 = f_q + f5 + f1
f10 = f_q + f6 + f2
f11 = f_q + f5 + f2
f12 = f_q + f6 + f1

f = [ f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 ]
# f = [ f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 ]
# f = [ f1 f2 f3 f4 f5 f6 f7 f8 ]
# f = [ f1 f2 f3 f4 f5 f6]
# f = [ f1 f2 f3 f4 ]
# f = f1
