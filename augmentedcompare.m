%% Longitudinal Augmented %%
s = tf('s');
kp = 1;
kd = 1.5552;
tauf = 20;
k = 1*pi/180;
num = (-0.3857*s^3 - 0.2926*s^2 + 0.000393*s - 1.487e-19);
den =  (s + 1.36 +4.33i)*(s + 1.36 -4.33i)*(s + 2.25e-2 - 9.64e-2i)*(s + 2.25e-2 + 9.64e-2i);
sys = ((kp + kd*s)/tauf*s + 1)*(num/den);
damp(sys*k)
