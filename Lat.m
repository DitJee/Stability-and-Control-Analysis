clear
close all
clc

%% Dimensionless Lateral-directional Body Axes State Space Equation

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
%Dimensionless to Body Axis Lateral A B
gam = 2;
alp = 8.4; 
Vo = 228.3151;
m = 23200; 
Ix = 36424;
Iz =  186248;
Ixz = 3600;
b = 12.787; 
rho = 0.8544; 
S =  54.148; 
g = 9.81; 
cbar =  54.148/12.787; 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mp = m/(0.5*rho*Vo*S);
Ixp = Ix/(0.5*rho*Vo*S*b);
Izp = Iz/(0.5*rho*Vo*S*b);
Ixzp = Ixz/(0.5*rho*Vo*S*b);
Ue = Vo*cos(alp*pi/180);
We = Vo*sin(alp*pi/180);
zta = gam + alp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = [mp 0 0 0 0
0 Ixp -Ixzp 0 0
0 -Ixzp Izp 0 0
0 0 0 1 0
0 0 0 0 1]
Ap = [Yv (Yp*b + mp*We) (Yr*b-mp*Ue) mp*g*cos(zta*pi/180) mp*g*sin(zta*pi/180)
Lv Lp*b Lr*b 0 0
Nv Np*b Nr*b 0 0
0 1 0 0 0
0 0 1 0 0]
Bp = [Vo*Yai Vo*Yrud
Vo*Lai Vo*Lrud
Vo*Nai Vo*Nrud
0 0
0 0]
A = M\Ap
B = M\Bp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transforming from Body Axes Dimensionless to Body Axes Dimensional
c = cbar; 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dimensionnal Lateral-directional Wind Axes State Space Equation
%Dimensional to Wind Axis Lateral A B
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ue = Vo;
We = 0;
zta = gam + alp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set transfer function from Numerator and denomenator
% input kk: 1 for aileron input
% input kk: 2 for aileron input
kk = 1;

[num,den]=ss2tf(A,B,C,D,kk);
[num1,den1]=ss2tf(A,B,C,D,1);
[num2,den2]=ss2tf(A,B,C,D,2);
% input step input
k = 1;
% if input is in degree -> input 0
% if input is in radian -> input 1
angle = 0
if(angle == 0)
    k = (pi/180)*k;
end
for i=1:size(C,1)
    G(i) = tf((num(i,:)),den);
end

%transfer function

for i=1:size(C,1)
    fprintf('\n For output %g : \n',i);
    G(i)
    fprintf('\n Factorised numerator');
    factor(poly2sym(num(i,:),sym('s')), 'FactorMode', 'real')
end
fprintf('\n Factorised denominator');
factor(poly2sym(den,sym('s')), 'FactorMode', 'real')
   
%step response plot
% input the desired variable (1 2 3 4 5)
stepresp = 1;
step(G(stepresp)*k)
        
%Final value

syms s
for i=1:size(C,1)
    trf = poly2sym(num(i,:),sym('s'))/poly2sym(den,sym('s'));
    fv = limit(trf*s*(k/s));
    if i==1
        fprintf('\n Variable %g : %.5g m/s',i,vpa(fv));
    else if i ==2
        fprintf('\n Variable %g : %.5g deg/s',i,vpa(fv/(pi/180)));
    else if i ==3
        fprintf('\n Variable %g : %.5g deg/s',i,vpa(fv/(pi/180)));
    else if i >3
        fprintf('\n Variable %g : %.5g deg',i,vpa(fv/(pi/180)));
        end
        end
        end
    end
end
   
%Eigen value and Eigen vector

[V,M] = eig(A);
EigenValue = M
fprintf('\nNOTES :');
fprintf('\n lamdas 0 0 0 ');
fprintf('\n 0 lamdas* 0 0 ');
fprintf('\n 0 0 lamdap 0 ');
fprintf('\n 0 0 0 lamdas*\n');
EigenVectorMagnitude = abs(V)
fprintf('\n Dutch Roll        Roll     Spiral');
    
%Damping Ratio and Undamped Natural Frequency
    
dd = factor(poly2sym(den,sym('s')), 'FactorMode', 'real');
d1 = sym2poly(dd(1));
d2 = sym2poly(dd(2));
d3 = sym2poly(dd(3));
if size(d1)==3
    d4 =d3;
    d3 =d1;
    d1 =d3;
else if size(d2) == 3
    d4 =d3;
    d3 =d2;
    d2 =d3;
    end
end
T1 = abs(1/roots(d1));
T2 = abs(1/roots(d2));
if T1 > T2
    Ts = -1/roots(d1);
    Tr = -1/roots(d2);
else
    Ts = -1/roots(d2);
    Tr = -1/roots(d1);
end
wd = sqrt(abs(d3(size(d3,2))));
drd = d3(size(d3,2)-1)/(2*wd);
fprintf('\n Tr = %.3g s',Tr);
fprintf('\n Ts = %.3g s\n',Ts);
fprintf('\n Oscillatory Dutch Roll mode:');
fprintf('\n Damping Ratio = %.3g',drd);
fprintf('\n Undamped Natural Frequency = %.3g rad/s\n',wd);
    
%SISO tool
% input the desired variable (1 2 3 4 5)
%siso = 1;
%sisotool(G(siso)*k)
