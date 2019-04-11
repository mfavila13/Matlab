%% function [Y] = Istim_rk(t,t1,ts, stim, on1)
% BME 445 2018.02.20
% This function creates the Action Potential using the Runge-Kutta 
% Method.  The input terms are:  
% t = initial time; t1 = final time; ts = timestep; stim = first   
% stimulation amplitude; on1 = switch to activate image saving 
% feature at the end of program.

function [Y] = Istim_rk(t,t1,ts, stim, on1)
%% Declaration of variables

% Potentials in units of mV
v = 0;                      % Resting Potential
e_K = -12.26;               
e_Na = 117.56;

% Conductivity in units of mS/cm^2
g_K = 36;
g_Na = 120;
g_L = 0.3;

% Number of entries
num = (t1-t)/ts;

% Modification to accommodate two stimuli in Part C
count = 0;
while count < 1 | count > 2
    count = input('How many stimuli?:  ');
    delay = input('What is delay of first stimulus?:  ');
    duration = input('What is duration of first stimulus?:  ');
    dur = duration/(t1-t)*num;
    D1 = delay/(t1-t)*num;
    I_stim = [zeros(1,D1) ones(1,dur)];
    I_stim = stim*[I_stim, zeros(1,num-length(I_stim)+1)];

    if count == 2
        x = input('Duration between Stimuli:  ');
        I_stim2 = [zeros(1,D1) zeros(1,dur) zeros(1,x/(t1-t)*...
        num) ones(1,dur)];
     
        % For Part 9
        stim = input('stim :  ');  % Ten times stimulus current.
        I_stim2 = 10*stim*[I_stim2, zeros(1,num-length(I_stim2)+1)];
        I_stim3 = I_stim + I_stim2;
        I_stim = I_stim3;
    end
end
        
% Calculation of initial values
m = alpha_m(v)/(alpha_m(v)+beta_m(v));
n = alpha_n(v)/(alpha_n(v)+beta_n(v));
h = alpha_h(v)/(alpha_h(v)+beta_h(v));

% Solving for e_L
e_L = v + (ik(v,n,g_K,e_K) + ina(v,m,h,g_Na,e_Na))/g_L;
fprintf('e_L was determined to be %2.3f mV at rest',e_L);
% e.g. e_L was determined to be 10.826 mV at rest

% Declaration of Variables:
T = t; V = v; 
M = m;  N = n; H = h; G_k = g_K; G_na = g_Na;
I_na = ina(v,m,h,g_Na,e_Na);
I_k = ik(v,n,g_K,e_K);
I_l = il(v,e_L);


% Iteration Process: Deriving equations using Runge-Kutta Method
for i = 1:num
% Runge-Kutta Algorithm for m
    m1 = f_m(m,v);
    m2 = f_m(m+ts*m1/2,v);
    m3 = f_m(m+ts*m2/2,v);
    m4 = f_m(m+ts*m3,v);
    m0 = (m1 + 2*m2 + 2*m3 + m4)/6;
    
% Runge-Kutta Algorithm for n
    n1 = f_n(n,v);
    n2 = f_n(n+ts*n1/2,v);
    n3 = f_n(n+ts*n2/2,v);
    n4 = f_n(n+ts*n3,v);
    n0 = (n1 + 2*n2 + 2*n3 + n4)/6;
    
% Runge-Kutta Algorithm for h
    h1 = f_h(h,v);
    h2 = f_h(h+ts*h1/2,v);
    h3 = f_h(h+ts*h2/2,v);
    h4 = f_h(h+ts*h3,v);
    h0 = (h1 + 2*h2 + 2*h3 + h4)/6;

% Solving for next timestep for m, n, h.
    m = m + ts*m0;
    n = n + ts*n0;
    h = h + ts*h0;

% Saving current m,n,h values.
    M = [M; m];
    N = [N; n];
    H = [H; h];


% Runge-Kutta Algorithm for V_m
    v1 = f_v(v,m,n,h,I_stim(i),e_L);
    v2 = f_v(v+ts*v1/2,m,n,h,I_stim(i),e_L);
    v3 = f_v(v+ts*v2/2,m,n,h,I_stim(i),e_L);
    v4 = f_v(v+ts*v3,m,n,h,I_stim(i),e_L);
    v0 = (v1 + 2*v2 + 2*v3 + v4)/6;
    v = v+ts*v0;
    t = t+ts;
    T = [T; t];
    V = [V; v];

    
% Saving current and conductance values 
    I_na = [I_na; ina(v,m,h,g_Na,e_Na)];
    G_na = [G_na; ina(v,m,h,g_Na,e_Na)/(v-e_Na)];
    
    I_k = [I_k; ik(v,n,g_K,e_K)];
    G_k = [G_k; ik(v,n,g_K,e_K)/(v-e_K)];    
end
%%%%%%%%%%%%%%%%
figure('NumberTitle', 'off', 'Name', 'Istim_rk Plots')
% Create Necessary Plots
subplot(4,2,1);
plot(T,V);
title('Membrane potential vs time plot of Action Potential ');
xlabel('time (ms)');
ylabel('Membrane Potential (mV)');
             
e_Na = 117.56;
axis([0,t1,min(e_K, e_Na)-10, max(e_K, e_Na)+10]);

subplot(4,2,2);
plot(T,-I_stim);
title('Stimulus current vs time plot of Action Potential');
xlabel('time (ms)');
ylabel('Stimulus current (uA)');

axis([0,t1,0,max(100,2*max(-I_stim(:)))]);


if on1 ~= 1
%{
pause
saveas(gcf, 'Part 3', 'jpg')

clf
%}
subplot(4,2,3);
plot(T,M);
title('m vs time plot of Action Potential ');
xlabel('time (ms)');
% ylabel('m value');
% subplot(3,1,2);
hold on
plot(T,H);
% title('h vs time plot of Action Potential');
% xlabel('time (ms)');
% ylabel('h value');
% subplot(3,1,3);
plot(T,N);
% title('n vs time plot of Action Potential');
% xlabel('time (ms)');
% ylabel('n value');
legend('m', 'h', 'n');
hold off
%{
pause
saveas(gcf, 'Part 4', 'jpg')

clf
%}
subplot(4,2,5);
plot(T,G_na);
title('G_N_a vs time plot of Action Potential ');
xlabel('time (ms)');
ylabel('Conductivity (mS/cm^2)');
subplot(4,2,6);
plot(T,G_k);
title('G_K vs time plot of Action Potential');
xlabel('time (ms)');
ylabel('Conductivity (mS/cm^2)');
%{
pause
saveas(gcf, 'Part 5', 'jpg')

clf
%}
subplot(4,2,7);
plot(T,I_na);
title('I_N_a vs time plot of Action Potential ');
xlabel('time (ms)');
ylabel('Current (uA/cm^2)');
subplot(4,2,8);
plot(T,I_k);
title('I_K vs time plot of Action Potential');
xlabel('time (ms)');
ylabel('Current (uA/cm^2)');
%{
pause
saveas(gcf, 'Part 6', 'jpg')
%}

end



%%%%%%%%%%%%%%
Y = [M,N,H,V,G_k,G_na,I_k,I_na];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f_v.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vp = f_v(v,m,n,h,I_stim, e_L)
% This is the differential equation for the membrane potential in 
% terms of the four currents and specific membrane capacitance C_m.

% Setting values into variables
C_m = 1; g_Na = 120;  g_K = 36;  e_K = -12.26;  e_Na = 117.56;

% Differential equation
vp = -1/C_m*(ina(v,m,h,g_Na,e_Na) + ik(v,n,g_K,e_K) + ...
il(v,e_L)+ I_stim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha_m.m, alpha_h.m, alpha_n.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpham = alpha_m(v_m)
% Alpha equation of m.

alpham = 0.1*(25-v_m)/(exp((25-v_m)/10)-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alphah = alpha_h(v_m)
% Alpha equation of h.

alphah = 0.07*exp(-v_m/20);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alphan = alpha_n(v_m)
% Alpha equation of n.

alphan = 0.01*(10-v_m)/(exp((10-v_m)/10)-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% beta_m.m, beta_h.m, beta_n.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  betam = beta_m(v_m)
% Beta equation of m

betam = 4.0*exp(-v_m/18.0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  betah = beta_h(v_m)
% Beta equation of h

betah = 1/(exp((30-v_m)/10)+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  betan = beta_n(v_m)
% Beta equation of n

betan = 0.125*exp(-v_m/80.0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f_m.m, f_h.m, f_n.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mp = f_m(m,v_m)
% Differential Equation of m

mp = alpha_m(v_m)*(1-m) - beta_m(v_m)*m;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hp = f_h(h, v_m)
% Differential Equation of h

hp = alpha_h(v_m)*(1-h) - beta_h(v_m)*h;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function np = f_n(n, v_m)
% Differential Equation of n

np = alpha_n(v_m)*(1-n) - beta_n(v_m)*n;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I_K = ik(v,n,g_K,e_K)

g_K = 36;
e_K = -12.26;
I_K = g_K*n.^4.*(v-e_K);
end

function I_Na = ina(v,m,h,g_Na,e_Na)

g_Na = 120;
e_Na = 117.56;
I_Na = g_Na.*m.^3.*h.*(v-e_Na);
end

function I_L = il(v, e_L)

% Declare variables 
g_L = 0.3;
I_L = g_L.*(v-e_L);
end