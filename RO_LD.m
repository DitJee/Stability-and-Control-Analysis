clear
close all
clc
%%Dimensionless Longitudinal Body Axes State Space Equation
%Dimensionless to Body Axis Longitudinal A B
gam = 2; %Flight Path Angle (deg)
alp = 8.4; %Body incidence (deg)
Vo = 228.3151; %m/s
m = 23200; %kg
Iy = 176790;
rho = 0.8544; %air density (kg/m^3)
S =  54.148; %Wind Area (m^3)
cbar =  54.148/12.787; %Mean ARD chord (m)
g = 9.81; %m/s^2
% de
Xu =  0.0072;
Xw =  0.0488;
Xwd = 0;
Xq = 0;
Xele =  0.0494;
Zu = -0.7231;
Zw = -3.8215;
Zwd = -0.5201;
Zq = 1.4291;
Zele = -0.4282;
Mu =  0.0612;
Mw = -0.2422;
Mwd = -0.6325;
Mq = -1.2458;
Mele = -0.5842;
%------------------------------------------------------
mp = m/(0.5*rho*Vo*S);
Iyp = Iy/(0.5*rho*Vo*S*cbar);
Ue = Vo*cos(alp*pi/180);
We = Vo*sin(alp*pi/180);
zta = gam + alp;
%-----------------------------------------------------
%------------------------------------------------------
M = [mp -Xwd*cbar/Vo 0 0
0 mp-(Zwd*cbar/Vo) 0 0
0 -Mwd*cbar/Vo Iyp 0
0 0 0 1]
Ap = [Xu Xw (Xq*cbar - mp*We) -mp*g*cos(zta*pi/180)
Zu Zw (Zq*cbar + mp*Ue) -mp*g*sin(zta*pi/180)
Mu Mw Mq*cbar 0
0 0 1 0]
Bp = [Vo*Xele
Vo*Zele
Vo*Mele
0]
A = M\Ap
B = M\Bp
%------------------------------------------------------
%% Dimensionless Lateral-directional Body Axes State Space Equation
%Dimensionless to Body Axis Lateral A B
gam = 2; %Flight Path Angle (deg)
alp = 8.4; %Body incidence (deg)
Vo = 228.3151; %m/s
m = 23200; %kg
Ix = 36424;
Iz =  186248;
Ixz = 3600;
b = 12.787; %Wing Span (m)
rho = 0.8544; %air density (kg/m^3)
S =  54.148; %Wind Area (m^3)
g = 9.81; %m/s^2
Yv = -0.6885;
Yp = 0;
Yr = 0;
Yai = -0.0214;
Yrud =  0.1452;
Lv = -0.1852;
Lp = -0.0958;
Lr =  0.0546;
Lai = 0.03256;
Lrud =  0.0082;
Nv =  0.0850;
Np = -0.0056;
Nr = -0.1250;
Nai =  0.0009;
Nrud = -0.0891;
%------------------------------------------------------
mp = m/(0.5*rho*Vo*S);
Ixp = Ix/(0.5*rho*Vo*S*b);
Izp = Iz/(0.5*rho*Vo*S*b);
Ixzp = Ixz/(0.5*rho*Vo*S*b);
Ue = Vo*cos(alp*pi/180);
We = Vo*sin(alp*pi/180);
zta = gam + alp;
%------------------------------------------------------
%-----------------------------------------------------
M = [mp 0 0 0 0
0 Ixp -Ixzp 0 0
0 -Ixzp Izp 0 0
0 0 0 1 0
0 0 0 0 1]
Ap = [Yv (Yp*b + mp*We) (Yr*b-mp*Ue) mp*g*cos(zta*pi/180) mp*g*sin(zta*pi/180)
Lv Lp*b Lr*b 0 0
Nv Np*b Nr*b 0 0
0 1 0 0 0
0 0 1 0 0];
Bp = [Vo*Yai Vo*Yrud
Vo*Lai Vo*Lrud
Vo*Nai Vo*Nrud
0 0
0 0]
A = M\Ap
B = M\Bp
%-----------------------------------------------------
%% Transforming from Body Axes Dimensionless to Body Axes Dimensional
%Body Axes
%Dimensionless to Dimensional
% gam = 0; %Flight Path Angle (deg)
% alp = 20; %Body incidence (deg)
% Vo = 263.2038583; %m/s
% m = 16500; %kg
% Iy = 172500;
% Ix = 42350;
% Iz = 187250;
% Ixz = 3250;
% rho = 0.54895; %air density (kg/m^3)
% S = 38.8; %Wind Area (m^3)
% b = 11.43; %Wing Span (m)
c = cbar; %Mean ARD chord (m)
% g = 9.81; %m/s^2
% Xu = 0.075;
% Xw = 0.0351;
% Xwd = 0;
% Xq = 0;
% Xele = 0.0627;
% Zu = -0.933;
% Zw = -2.15;
% Zwd = -0.2206;
% Zq = 0.063;
% Zele = -0.3747;
% Mu = 0.0421;
% Mw = -0.2151;
% Mwd = -0.5630;
% Mq = -1.4452;
% Mele = -0.5445;
% Yv = -0.5859;
% Yp = 0;
% Yr = 0;
% Yai = -0.0211;
% Yrud = 0.1814;
% Lv = -0.1214;
% Lp = -0.1142;
% Lr = 0.04201;
% Lai = 0.0730;
% Lrud = 0.0062;
% Nv = 0.0411;
% Np = -0.0052;
% Nr = -0.0911;
% Nai = 0.0008;
% Nrud = -0.0829;
fprintf('Dimensional Value:...\n');
DXu = Xu*0.5*rho*Vo*S
DXw = Xw*0.5*rho*Vo*S
DXwd = Xwd*0.5*rho*S*c
DXq = Xq*0.5*rho*Vo*S*c
DXele = Xele*0.5*rho*(Vo^2)*S
DZu = Zu*0.5*rho*Vo*S
DZw = Zw*0.5*rho*Vo*S
DZwd = Zwd*0.5*rho*S*c
DZq = Zq*0.5*rho*Vo*S*c
DZele = Zele*0.5*rho*(Vo^2)*S
DMu = Mu*0.5*rho*Vo*S*c
DMw = Mw*0.5*rho*Vo*S*c
DMwd = Mwd*0.5*rho*S*(c^2)
DMq = Mq*0.5*rho*Vo*S*(c^2)
DMele = Mele*0.5*rho*(Vo^2)*S*c
DYv = Yv*0.5*rho*Vo*S
DYp = Yp*0.5*rho*Vo*S*b
DYr = Yr*0.5*rho*Vo*S*b
DYai = Yai*0.5*rho*(Vo^2)*S
DYrud = Yrud*0.5*rho*(Vo^2)*S
DLv = Lv*0.5*rho*Vo*S*b
DLp = Lp*0.5*rho*Vo*S*b*b
DLr = Lr*0.5*rho*Vo*S*b*b
DLai = Lai*0.5*rho*(Vo^2)*S*b
DLrud = Lrud*0.5*rho*(Vo^2)*S*b
DNv = Nv*0.5*rho*Vo*S*b
DNp = Np*0.5*rho*Vo*S*b*b
DNr = Nr*0.5*rho*Vo*S*b*b
DNai = Nai*0.5*rho*(Vo^2)*S*b
DNrud = Nrud*0.5*rho*(Vo^2)*S*b
%-----------------------------------------------------
%%  Transforming Body Axes Dimensional to Wind Axes Dimensional
%Dimensional Body to Wind Axes
co = cos(8.4*pi/180);
co2 = (cos(8.4*pi/180))^2;
sn = sin(8.4*pi/180);
sn2 = (sin(8.4*pi/180))^2;
sc = sn*co;
WDXu = DXu*co2 + DZw*sn2 + (DXw + DZu)*sc
WDXw = DXw*co2 - DZu*sn2 - (DXu - DZw)*sc
WDXwd = DXwd*co2 + DZwd*sc
WDXq = DXq*co + DZq*sn
WDXele = DXele*co + DZele*sn
WDYv = DYv
WDYp = DYp*co + DYr*sn
WDYr = DYr*co - DYp*sn
WDYai = DYai
WDYrud = DYrud
WDZu = DZu*co2 - DXw*sn2 - (DXu - DZw)*sc
WDZw = DZw*co2 + DXu*sn2 - (DXw + DZu)*sc
WDZwd = DZwd*co2 - DXwd*sc
WDZq = DZq*co - DXq*sn
WDZele = DZele*co - DXele*sn
WDLv = DLv*co + DNv*sn
WDLp = DLp*co2 + DNr*sn2 + (DLr + DNp)*sc
WDLr = DLr*co2 - DNp*sn2 - (DLp - DNr)*sc
WDLai = DLai*co + DNai*sn
WDLrud = DLrud*co + DNrud*sn
WDMu = DMu*co + DMw*sn
WDMw = DMw*co - DMu*sn
WDMwd = DMwd*co
WDMq = DMq
WDMele = DMele
WDNv = DNv*co - DLv*sn
WDNp = DNp*co2 - DLr*sn2 - (DLp - DNr)*sc
WDNr = DNr*co2 + DLp*sn2 - (DLr + DNp)*sc
WDNai = DNai*co - DLai*sn
WDNrud = DNrud*co - DLrud*sn
%------------------------------------------------------
%% Dimensionnal Lateral-directional Wind Axes State Space Equation
%Dimensional to Wind Axis Lateral A B
% gam = 0; %Flight Path Angle (deg)
% alp = 20; %Body incidence (deg)
% Vo = 263.2038583; %m/s
% m = 16500; %kg
% Ix = 42350;
% Iz = 187250;
% Ixz = 3250;
% rho = 0.3809; %air density (kg/m^3)
% S = 49.239; %Wing Area (m^3)
% b = 11.787; %Wing Span (m)
% g = 9.81; %m/s^2
Yv = WDYv
Yp = WDYp
Yr = WDYr
Yai = WDYai
Yrud = WDYrud
Lv = WDLv
Lp = WDLp
Lr = WDLr
Lai = WDLai
Lrud = WDLrud
Nv = WDNv
Np = WDNp
Nr = WDNr
Nai = WDNai
Nrud = WDNrud
%------------------------------------------------------
Ue = Vo;
We = 0;
zta = gam + alp;
%------------------------------------------------------
M = [ m 0 0 0
    0 Ix -Ixz 0
0 -Ixz Iz 0
0 0 0 1]
Ap = [Yv Yp (Yr-m*Ue) m*g
Lv Lp Lr 0
Nv Np Nr 0
0 1 0 0]
Bp = [Yai Yrud
Lai Lrud
Nai Nrud
0 0]
A = M\Ap
B = M\Bp
%------------------------------------------------------
%%  Full-model of Lateral-directional Analysis
%Setting matrix A,B,C,D
%lateral-directional transfer function
%STATE EQUATION MATRIX : xdot = Ax +Bu(input)
A=A;
B=B;

%OUTPUT EQUATION : y = Cx + Du
C= [1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1
1/Vo 0 0 0]
D=[0 0;0 0;0 0;0 0;0 0]

%-------------------------------------------------------------------------%
%set transfer function from Numerator and denomenator
kk = input('\nDetermine for the aileron input(1) or rudder input(2)....');
[num,den]=ss2tf(A,B,C,D,kk);
[num1,den1]=ss2tf(A,B,C,D,1);
[num2,den2]=ss2tf(A,B,C,D,2);
k = input('\nStep input = ');
angle = input('\nDEGREE(0) or RAD(1) = ');
if(angle == 0)
    k = (pi/180)*k;
end
for i=1:size(C,1)
    G(i) = tf((num(i,:)),den);
end
%% reduce shit
%Setting matrix A,B,C,D
%STATE EQUATION MATRIX : xdot = Ax +Bu(input)
%CONCISE VALUE
yv = -0.156734375998225
yp = 0
yr = -2.283151000000000e+02
yphi = 9.81
lv = -0.313285807442403
lp = -2.127174695208387
lr = -1.111511594721804
lphi = 0
nv = 0.034244602292882
np = -0.091494796349657
nr = -0.588026461515156
nphi = 0
yai = -1.112266157168786
yrud= 7.546777851444290
lai = 13.685451407422224
lrud= 2.812318440543014
nai = -0.055535740469372
nrud= -7.450682590007344
Ar=[lp 0
1 0]
Br=[lai lrud
0 0]
%OUTPUT EQUATION : y = Cx + Du
Cr=[1 0
0 1 ];
Dr=[0 0;0 0];
%-------------------------------------------------------------------------%
%set transfer function from Numerator and denomenator
[num,den]=ss2tf(Ar,Br,Cr,Dr,1);
k = input('\nStep input = ');
angle = input('\nDEGREE(0) or RAD(1) = ');
if(angle == 0)
    k = 0.01745*k;
end
for i=1:size(Cr,1)
    Gr(i) = tf((num(i,:)),den);
end
cond = 0;
while cond == 0
    fprintf('\nWHICH VALUE DO YOU WANT TO FIND');
    fprintf('\n1) Transfer Function of roll mode reduced order ');
    fprintf('\n2) Step Response Plot of roll mode reduced order');
    fprintf('\n3) Approximated Characteristics of three modes');
    choice = input('\n.....');
    %Transfer Function of roll mode reduced order
    if choice == 1
        for i=1:size(Cr,1)
            fprintf('\n For output %g : \n',i);
            Gr(i)
            fprintf('\n Factorised numerator');
            factor(poly2sym(num(i,:),sym('s')), 'FactorMode', 'real')
        end
        fprintf('\n Factorised denominator');
        factor(poly2sym(den,sym('s')), 'FactorMode', 'real')
    end
    %Step Response Plot of roll mode reduced order
    if choice == 2
        cont = 0;
        while(cont == 0)
            stepresp = input('\n Step response of variable....');
            if stepresp == 1
                step(Gr(1),5);
                hold on
                step(G(2),5);
                hold off
                title('$p(s)/\xi(s)$', 'Interpreter', 'latex', 'FontSize', 16);
                legend("Reduced order Model, Roll Mode","Full-order Model",'Location','southeast')
            end
            if stepresp == 2
                step(Gr(2),5);
                hold on
                step(G(4),5);
                hold off
                title('$\phi(s)/\xi(s)$', 'Interpreter', 'latex', 'FontSize', 16);
                legend("Reduced order Model, Roll Mode","Full-order Model",'Location','southeast')
            end
            cont = input('\n Do you want to find more step response? YES(0) or NO(1)');
        end
    end
    %Approximated Characteristics of three modes
    if choice == 3
        dd = factor(poly2sym(den,sym('s')), 'FactorMode', 'real');
        d1 = sym2poly(dd(1));
        d2 = sym2poly(dd(2));
        if size(d1) > size(d2)
            Tr = abs(1/roots(d1));
        else
            Tr = abs(1/roots(d2));
        end
        Ts = yr*((lv*np)-(lp*nv))/(yphi*((lr*nv)-(lv*nr)));
        wd = sqrt(nr*yv-nv*yr);
        drd = -(nr+yv)/(2*wd);
        fprintf('\n Tr (Roll Subsidence Time Constant) = %.5g s\n',Tr);
        fprintf('\n Ts (Spiral Time Constant) = %.5g s\n',Ts);
        fprintf('\n Oscillatory Dutch Roll mode:');
        fprintf('\n Damping Ratio = %.5g',drd);
        fprintf('\n Undamped Natural Frequency = %.5g rad/s\n',wd);
    end
end
