clear
clc
gam = 2; 
alp = 8.4;
Mach = 0.7;
Ta = 264.772;
V0 = sqrt(1.4*287*Ta)
rho = 0.8546;
m = 23200; 
Iy = 176790;
rho = 0.54895; 
S = 54.148; 
b = 12.787
c = S/b; 
Ix = 36424;
Iz = 186248;
Ixz = 3250;
g = 9.81; 
Xu = 0.0072;
Xw = 0.0488;
Xwd = 0;
Xq = 0;
Xele = 0.0494;
Zu = -0.7231;
Zw = -3.8215;
Zwd = -0.5201;
Zq = -1.4291;
Zele = -0.4282;
Mu = 0.0612;
Mw = -0.2422;
Mwd = -0.6325;
Mq = -1.2458;
Mele = -0.5842;
Yv = -0.6885;
Yp = 0;
Yr = 0;
Yai = -0.0214;
Yrud = 0.1452;
Lv = -0.1852;
Lp = -0.0958;
Lr = 0.0546;
Lai = 0.03256;
Lrud = 0.0082;
Nv = 0.0850;
Np = -0.0056;
Nr = -0.1250;
Nai = 0.0009;
Nrud = -0.0891;
DXu = Xu*0.5*rho*V0*S;
DXw = Xw*0.5*rho*V0*S;
DXwd = Xwd*0.5*rho*S*c;
DXq = Xq*0.5*rho*V0*S*c;
DXele = Xele*0.5*rho*(V0^2)*S;
DZu = Zu*0.5*rho*V0*S;
DZw = Zw*0.5*rho*V0*S;
DZwd = Zwd*0.5*rho*S*c;
DZq = Zq*0.5*rho*V0*S*c;
DZele = Zele*0.5*rho*(V0^2)*S;
DMu = Mu*0.5*rho*V0*S*c;
DMw = Mw*0.5*rho*V0*S*c;
DMwd = Mwd*0.5*rho*S*(c^2);
DMq = Mq*0.5*rho*V0*S*(c^2);
DMele = Mele*0.5*rho*(V0^2)*S*c;
DYv = Yv*0.5*rho*V0*S;
DYp = Yp*0.5*rho*V0*S*b;
DYr = Yr*0.5*rho*V0*S*b;
DYai = Yai*0.5*rho*(V0^2)*S;
DYrud = Yrud*0.5*rho*(V0^2)*S;
DLv = Lv*0.5*rho*V0*S*b;
DLp = Lp*0.5*rho*V0*S*b*b;
DLr = Lr*0.5*rho*V0*S*b*b;
DLai = Lai*0.5*rho*(V0^2)*S*b;
DLrud = Lrud*0.5*rho*(V0^2)*S*b;
DNv = Nv*0.5*rho*V0*S*b;
DNp = Np*0.5*rho*V0*S*b*b;
DNr = Nr*0.5*rho*V0*S*b*b;
DNai = Nai*0.5*rho*(V0^2)*S*b;
DNrud = Nrud*0.5*rho*(V0^2)*S*b;
cost = 0.93;
cost2 = 0.88;
sint = 0.34;
sint2 = 0.121;
sc = sint*cost;
WDXu = DXu*cost2 + DZw*sint2 + (DXw + DZu)*sc;
WDXw = DXw*cost2 - DZu*sint2 - (DXu - DZw)*sc;
WDXwd = DXwd*cost2 + DZwd*sc;
WDXq = DXq*cost + DZq*sint;
WDXele = DXele*cost + DZele*sint;
WDYv = DYv;
WDYp = DYp*cost + DYr*sint;
WDYr = DYr*cost - DYp*sint;
WDYai = DYai;
WDYrud = DYrud;
WDZu = DZu*cost2 - DXw*sint2 - (DXu - DZw)*sc;
WDZw = DZw*cost2 + DXu*sint2 - (DXw + DZu)*sc;
WDZwd = DZwd*cost2 - DXwd*sc;
WDZq = DZq*cost - DXq*sint;
WDZele = DZele*cost - DXele*sint;
WDLv = DLv*cost + DNv*sint;
WDLp = DLp*cost2 + DNr*sint2 + (DLr + DNp)*sc;
WDLr = DLr*cost2 - DNp*sint2 - (DLp - DNr)*sc;
WDLai = DLai*cost + DNai*sint;
WDLrud = DLrud*cost + DNrud*sint;
WDMu = DMu*cost + DMw*sint;
WDMw = DMw*cost - DMu*sint;
WDMwd = DMwd*cost;
WDMq = DMq;
WDMele = DMele;
WDNv = DNv*cost - DLv*sint;
WDNp = DNp*cost2 - DLr*sint2 - (DLp - DNr)*sc;
WDNr = DNr*cost2 + DLp*sint2 - (DLr + DNp)*sc;
WDNai = DNai*cost - DLai*sint;
WDNrud = DNrud*cost - DLrud*sint;
Xu = DXu;
Xw = DXw;
Xwd = DXwd;
Xq = DXq;
Xele = DXele;
Zu = DZu;
Zw = DZw;
Zwd = DZwd;
Zq = DZq;
Zele = DZele;
Mu = DMu;
Mw = DMw;
Mwd = DMwd;
Mq = DMq;
Mele = DMele;
Ue = V0*cos(alp*pi/180);
We = V0*sin(alp*pi/180);
zta = gam + alp;
M = [ m -Xwd 0 0
0 m-Zwd 0 0
0 -Mwd Iy 0
0 0 0 1];
Ap = [Xu Xw (Xq - m*We) -m*g*cos(zta*pi/180)
    Zu Zw (Zq + m*Ue) -m*g*sin(zta*pi/180)
Mu Mw Mq 0
0 0 1 0];
Bp = [Xele
Zele
Mele
0];
A = M\Ap;
B = M\Bp;
Xu = WDXu;
Xw = WDXw;
Xwd = WDXwd;
Xq = WDXq;
Xele = WDXele;
Zu = WDZu;
Zw = WDZw;
Zwd = WDZwd;
Zq = WDZq;
Zele = WDZele;
Mu = WDMu;
Mw = WDMw;
Mwd = WDMwd;
Mq = WDMq;
Mele = WDMele;
Ue = V0;
We = 0;
zta = gam + alp;
M = [ m -Xwd 0 0
    0 m-Zwd 0 0
    0 -Mwd Iy 0
    0 0 0 1];
Ap = [Xu Xw Xq -m*g;
      Zu Zw (Zq + m*Ue) 0;
      Mu Mw Mq 0;
       0 0 1 0];
Bp = [Xele
    Zele 
    Mele 
    0];
Ar = M\Ap;
Br = M\Bp;
% Longitudinal Analysis
s = tf('s');
A=A;
B=B;
C= [1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
0 1/263.2038583 0 0
0 -1/263.2038583 0 1];
D=[0;0;0;0;0;0];
[num,den]=ss2tf(A,B,C,D);
k = 1;
k = (pi/180)*k;

for i=1:size(C,1)
    G(i) = tf((num(i,:)),den);
end
for i=1:size(C,1)
    fprintf('\n For output %g : \n',i);
    G(i)
    fprintf('\n Factorised numerator');
    vpa(factor(poly2sym(num(i,:),sym('s')), 'FactorMode', 'real'),4)
end
fprintf('\n Factorised denominator');
factor(poly2sym(den,sym('s')), 'FactorMode', 'complex')
step(G(1)*k)
title('$u(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
hold on
figure()
step(G(2)*k)
title('$w(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
hold on
figure()
step(G(3)*k)
title('$q(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
hold on
figure()
step(G(4)*k)
title('$\theta(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
hold on
figure()
step(G(5)*k)
title('$\alpha(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
hold on
figure()
step(G(6)*k)
title('$\gamma(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
hold off
syms s
for i=1:size(C,1)
    trf = poly2sym(num(i,:),sym('s'))/poly2sym(den,sym('s'));
    fv = limit(trf*s*(k/s));
    if i<=2
        fprintf('\n Variable %g : %.8g m/s',i,vpa(fv));
    elseif i ==3
        fprintf('\n Variable %g : %.8g deg/s',i,vpa(fv/(pi/180)));
    else
        fprintf('\n Variable %g : %.8g deg',i,vpa(fv/(pi/180)));
    end
end
[V,M] = eig(A);
Eigen = M;
EigenMag = abs(V);
damp(G(1));;
sisotool(G(1)*k);
% G(i) to show each variable 
Ar= M\Ap;
Br = M\Bp;
Xu = Ar(1,1);
Xw = Ar(1,2);
Xq = Ar(1,3);
Xtheta = Ar(1,4);
Zu = Ar(2,1);
Zw = Ar(2,2);
Zq = Ar(2,3);
Ztheta = Ar(2,4);
Mu = Ar(3,1);
Mw = Ar(3,2);
Mq = Ar(3,3);
Mtheta = Ar(3,4);
Xele = Br(1);
Zele = Br(2);
Mele = Br(3);
%%%%% Reduced Order Short Period%%%%%
s = tf('s');
k = 1*pi/180;
wn = Zele*(s + Ue*(Mele/Zele))/(s^2 - (Mq + Zw)*s + (Mq*Zw - Mw*Ue));
wnreal = (-29.14*s^3 - 7119*s^2 - 111.7*s - 34.31)/(s^4 + 1.714*s^3 + 9.857*s^2 + 0.4558*s + 0.09582);
qn = Mele*(s - Zw)/(s^2 - (Mq + Zw)*s + (Mq*Zw - Mw*Ue));
qnreal = (-22.1*s^3 - 16.76*s^2 + 0.02252*s)/(s^4 + 1.714*s^3 + 9.857*s^2 + 0.4558*s + 0.09582);
step(wn*k,wnreal*k,8)
title('$w(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
legend('Reduced Order Model','Full Order Model')
hold on
figure()
step(qn*k,qnreal*k,8)
title('$q(s)/\eta(s)$', 'Interpreter', 'latex', 'FontSize', 16);
legend('Reduced Order Model','Full Order Model')
hold off
