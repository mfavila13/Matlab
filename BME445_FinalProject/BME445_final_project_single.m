%% BME 445 Lab Module 3
% This is a MATLAB Live Script (only available in advanced versions of MATLAB).  
%% Lab 5 Task 1 - study the overall code.  See insructions for details.
%% Set up Membrane Ring (George R. Mines Experiment)
% This part of the code sets up the membrane ring for a re-entry model.
%%
function BME445_final_project_single()
%num_x = 0.002*(1:30);
min_dx = 7;
max_dx = 30;
lag_pot = zeros(1,max_dx);
len_of_max = zeros(1,max_dx);
%for deltax = min_dx:max_dx
deltax = input('Input dx value (enter integer): ')
NumLocs = 100;                  % Number of locations
radius = 1.5/(2*pi);            % radius
d_theta = 2*pi./NumLocs;        % distance in radial coordinates
theta = d_theta*(0:NumLocs-1);
x=radius*cos(theta);            % Polar coordinates x
y=radius*sin(theta);            % Polar coordinates y
ground = zeros(1,NumLocs);      % Create the ring at rest potential
dx = deltax * 0.001

R_i = .5;
a = 0.002;

% Test out 3D plot of a ring
%plot3(x,y,ground,':b', 'LineWidth', 2);
%title('Membrane Potential around ring');

%% Set up Initial Variables and Input Options
% This is a simple interface using the input command to set up the 
%%
disp('BME 445 Final Project Tasks')
%{
mode_flag = input(...
    ['   Healthy mode         (press 1) ',...
    '\n   Unidirectional block (press 2) ',... 
    '\n   VFib model           (press 3) ', ...
    '\n   Cardiac Ablation sim    (press 4) ', ...
    '\n   Your choice: ',]);
%}
mode_flag = 1;
%'\n   _________________    (press 4) ]);
% Lab 6 task - Consider adding later
% Uncomment below for Lab #6 (Apr 3), and add to above options.
%'\n   VFib model           (press 3) ',... 
%'\n   _________________    (press 4) ',...

% Uncomment below for Lab #6 (Apr 3)
% add_defib = input('Add Defibrillator? \n  Yes (press 1), \n  No  (press 0): ');
add_defib = 0; %input('Add Defibrillator? \n  Yes (press 1), \n  No  (press 0): ');
dur_total = 20; %input('Total duration of simulation? (ms, deault 40) ');

dt = 0.01;
numsteps = dur_total/dt;
%% Set up regular stimulus currents for simulated heart beats
%%
% Create first Stimuli 
I_stim1 = 250;
tdur = 1;
tdelay = 0;
Nstim_pos1 = 75;

% Create second Stimuli
I_stim2 = 350;
tdur2 = 1;
tdelay2 = 20; %input('Time of 2nd stimuli: ');
Nstim_pos2 = 75; % input('Location of 2nd stimuli: ');

%% 
% Lab #6-7 - Debrillator code will be implemented in April.

%  This is the defibrillator variable.
I_defib = 350;
tdur3 = 1;
tdelay3 = dur_total;
%Nstim3=75;
%Nstim_pos3 = 100;

% Heart Stimulus after defibrillation that stabilizes the heart
I_stim4 = 350;
tdur4 = 1;
tdelay4 = 35;
Nstim4 = 1:100;

% Setting up the cable parameters.
v_m = ground;                      

C_m = 1; 
g_L_sup = 0.3; 
g_Na_sup = 120;  
g_K_sup = 36;  
e_K = -12.26;  
e_Na = 117.56;
%% Initialize m, n, h parameters.
%%
% Calculation of initial alpha and beta values.
m = alpha_m(v_m)./(alpha_m(v_m)+beta_m(v_m));
n = alpha_n(v_m)./(alpha_n(v_m)+beta_n(v_m));
h = alpha_h(v_m)./(alpha_h(v_m)+beta_h(v_m));
%% Initialize conductivity and current values.
%%
% Calculation of conductivity and current values

g_K = g_K_sup*n.^4;
g_Na = g_Na_sup*m.^3.*h;
g_L = g_L_sup*ones(1,NumLocs);

% Finding e_L
e_L = (ik(v_m,n,g_K,e_K) + ina(v_m,m,h,g_Na,e_Na))./g_L;

V = v_m; 
M = m;  N = n; H = h; G_k = g_K; G_na = g_Na;
I_na = ina(v_m,m,h,g_Na,e_Na);
I_k = ik(v_m,n,g_K,e_K);
I_l = il(v_m,g_L,e_L);

I_Na = I_na;
G_Na = I_na./(v_m-e_Na);

I_K = I_k;
G_K = G_k;

Iion = I_k + I_na + I_l;
I_S = ground;
%% Lab #5 Task 2. This part of the code sets up the Cable Equations. 
%%
% Matrix A will play a very important part of setting up 
% the cable equations for this Final Project.

% Setting up Matrix A:
b1 = 1; b2 = 2;
A = [-b2, b1, zeros(1,NumLocs-3), b1; b1, -b2, b1, zeros(1,NumLocs-3)];

for k = 3:NumLocs-1
    B = 0*linspace(1,NumLocs,NumLocs);
    B(k-1) = b1; B(k) = -b2; B(k+1) = b1;
    A = [A; B];
end

A = [A; [b1, [zeros(1,NumLocs-3)], b1,-b2]];
% Matrix A is set.

I_s = I_S;
coeff_A = a/(2*R_i*dx*dx)*ones(1,NumLocs);
%% Lab 5 Task 3 - Modification for Unidirectional Block (Disease Model)
%%
% unidirectional block (Disease)
if mode_flag == 2 || mode_flag == 3
    % In Lab #5, you are tasked to comment on the purpose of this modification.
    j = 10;
    A(j,j-1) = 0;
    A(j,j) = -1;
end
if mode_flag == 4
    sizeCA = 5;
    for i = 0:sizeCA
        j = 10;
        A(j+i,:) = 0;
    end
end
%% Implementing Forward Euler's Method
%%
%count = 1*dt;
%bool = 0;
pass = 0;
pass2 = 0;
node = 0;
time_until_0 = 0;
for k = 1:numsteps-1
% Solving for the present v_m
    v_0 = dt/C_m.*(coeff_A.*(A*v_m')' - Iion + I_s);
    v_m = v_m + v_0;
    V = [V; v_m];
    
% Runge-Kutta Algorithm for m, n, h
    m0 = rk_fcn('f_m',m,v_m,dt);
    n0 = rk_fcn('f_n',n,v_m,dt);
    h0 = rk_fcn('f_h',h,v_m,dt);
    
% Solving for next timestep for m, n, h.
    m = m + dt*m0;
    n = n + dt*n0;
    h = h + dt*h0;

    M = [M; m];
    N = [N; n];
    H = [H; h];

% solving for ion currents
    g_Na = g_Na_sup*m.^3.*h;
    I_na = g_Na.*(v_m-e_Na);
    I_Na = [I_Na; I_na];
%    G_Na = [G_Na; g_Na];

    g_K = g_K_sup*n.^4;
    I_k = g_K.*(v_m-e_K);
    I_K = [I_K; I_k];
%    G_K = [G_K; g_K];    

% Code for two stimuli
    
    if mode_flag == 1 || mode_flag == 2 ||  mode_flag == 3 || mode_flag ==4
    %  healthy (1)    or unidirectional block (2)
        I_s = 0*linspace(1,2,NumLocs);
        if k > tdelay/dt && k <= (tdelay+tdur)/dt
            I_s(Nstim_pos1) = I_stim1;
        %elseif k > tdelay2/dt && k <= (tdelay2+tdur2)/dt            
         %   I_s(Nstim_pos2) = I_stim2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Lab 6 (April 3rd):  Your task is to 
        % develop a dedicated model for 
        % application of defibrillator here.
    if add_defib == 1
        for o = 1:99
            if (V(k,o)) >= 50 && (V(k,o+1)) < -8 %&& k > 300
                %count
                %count = count + .01;
                %spot = o;
                if pass == 0
                    node = o + 1
                    time = k*.01;
                    tdelay3 = time;
                    pass = 1;
                end    
            end
        end
    end
    if k > tdelay3/dt && k <= (tdelay3+tdur3)/dt
        I_s(node) = I_defib;
        pass = 1;
    elseif pass ==1 && (tdelay3+tdur3)/dt + 1
        pass = 0;
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
        
    I_S = [I_S; I_s];
    Iion = I_k+I_na+I_l;
    waitbar(k/numsteps);
end
%% Simulating the ring model
%%

figure(445)
for j = 1:1
for k = 1:numsteps
    if mod(k,10) == 0
        P = plot3(x,y,V(k,:),x,y,ground,':b', 'LineWidth', 2);
        axis([-2*radius 2*radius -2*radius 2*radius -120 120]);
        title(sprintf('t  = %2.2f ms', k*dt));
        %pause(0.05)
        drawnow
    end
end
%V(500,:)
%size(V,2)
%size(V,1)
end
%This variable will be useful in obtaining extracellular potential measures

Time = linspace(0,dur_total-dt,numsteps); 
%% Obtaining measurements from a probe
% For Lab #7 (April 10-17), you will program an 'ECG Probe' to detect the measurements 
% during this AP propagation from a point source.  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

probe = [0,0,1];
j = 10;
nodes = (1:100);
for k = 1:numsteps
for i = 1:100
    Theta = 2*pi * i / 100; % position of 1 node in respect to 0
    X = radius * cos(Theta) ;
    Y = radius * sin(Theta);
    Z(k,i) = V(k,i) / 100; 
    distance(k,i) = ( (X - 0)^2 + (Y - 0)^2 + (Z(k,i) - 1)^2 )^ (1/2);
    flux_e(k,i) = distance(k,i)/ (0.01 + distance(k,i)) * (V(k,i)); 
    flux_i(k,i) = .01/ (0.01 + distance(k,i)) * (V(k,i)); 
    
end
%maxflux(k) = max(flux(k,:))
waitbar(k/numsteps)
end
ex_pot = flux_e - flux_i;
%sumflux = ( sum(flux') *100 );
maxpot = max(V');
pot_count = 0;
for i = 1:numsteps
    if maxpot(i) > 100
        pot_count = pot_count + 1;
    end
end
disp(pot_count * 0.01)
len_of_max(deltax) = pot_count * 0.01
maxcur_na = max(I_Na');
maxcur_k  = max(I_K');
timelength(deltax) = 0;
for i = 1:numsteps
    if maxpot(i) > 100
        timelength(deltax) = timelength(deltax) + 0.01;
    end
end
timelength(deltax)
[peak, peak_t] = max(ex_pot(:,90));
peak_t * 0.01
lag_pot(deltax) = peak_t * 0.01


figure(1)
subplot(4,1,1)
plot(Time,ex_pot(:,90))
title('Contractile Activity (at node 90)');...
    xlabel('Time (ms)');...
    ylabel('Potemtial (mV)');


subplot(4,1,2)
plot(Time,maxpot)
title('Maximum Potential (all nodes)');...
    xlabel('Time (ms)');...
    ylabel('Potemtial (mV)');

subplot(4,1,3)
plot(Time,maxcur_na)
title('Maximum Sodium Current');...
    xlabel('Time (ms)');...
    ylabel('Current (uA)');

subplot(4,1,4)
plot(Time,maxcur_k)
title('Maximum Potassium Current');...
    xlabel('Time (ms)');...
    ylabel('Current (uA)');
%figure('NumberTitle','off','Name','Membrane Potential at given node')
%plot(Time,(V(:,90)))



% Lab 7 (April 10th):  Your task is to 
% develop a point source to detect the 
% AP during this threshold stimuli.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%waitbar(deltax/30)
%figure('NumberTitle','off','Name','dx change')
%plot(num_x, timelength)
end
%{
figure(2)
subplot(2,1,1)
plot(min_dx:max_dx,len_of_max(min_dx:max_dx))
title('Time duration of Propogation');...
    xlabel('delta x value (*10^-4)');...
    ylabel('time (ms)');

subplot(2,1,2)
plot((min_dx:max_dx),lag_pot(min_dx:max_dx))
title('Membrane Potential Lag behind Stimulus Current');...
    xlabel('delta x value (*10^-4)');...
    ylabel('time (ms)');
%}
%end

%% 
% Runge kutta function
%%
function x0 = rk_fcn(fcn_str,x,v_m,dt)
    % m0 = rk_fcn('fcn_m',m,v_m,dt);
    x1 = eval(sprintf('%s(x,v_m);',fcn_str));
    x2 = eval(sprintf('%s(x+dt*x1/2,v_m);',fcn_str));
    x3 = eval(sprintf('%s(x+dt*x2/2,v_m);',fcn_str));
    x4 = eval(sprintf('%s(x+dt*x3,v_m);',fcn_str));
    x0 = (x1 + 2*x2 + 2*x3 + x4)/6;
end

% f_v.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vp = f_v(v,m,n,h,I_stim, e_L)
% This is the differential equation for the membrane potential in 
% terms of the four currents and specific membrane capacitance C_m.

% Setting values into variables
C_m = 1; g_Na = 120;  g_K = 36;  e_K = -12.26;  e_Na = 117.56;

% Differential equation
vp = -1./C_m*(ina(v,m,h,g_Na,e_Na) + ik(v,n,g_K,e_K) + ...
il(v,e_L)+ I_stim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alpha_m.m, alpha_h.m, alpha_n.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpham = alpha_m(v_m)
% Alpha equation of m.

alpham = 0.1.*(25-v_m)./(exp((25-v_m)./10)-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alphah = alpha_h(v_m)
% Alpha equation of h.

alphah = 0.07.*exp(-v_m./20);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alphan = alpha_n(v_m)
% Alpha equation of n.

alphan = 0.01.*(10-v_m)./(exp((10-v_m)./10)-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% beta_m.m, beta_h.m, beta_n.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  betam = beta_m(v_m)
% Beta equation of m

betam = 4.0.*exp(-v_m./18.0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  betah = beta_h(v_m)
% Beta equation of h

betah = 1./(exp((30-v_m)./10)+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  betan = beta_n(v_m)
% Beta equation of n

betan = 0.125.*exp(-v_m./80.0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f_m.m, f_h.m, f_n.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mp = f_m(m,v_m)
% Differential Equation of m

mp = alpha_m(v_m).*(1-m) - beta_m(v_m).*m;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hp = f_h(h, v_m)
% Differential Equation of h

hp = alpha_h(v_m).*(1-h) - beta_h(v_m).*h;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function np = f_n(n, v_m)
% Differential Equation of n

np = alpha_n(v_m).*(1-n) - beta_n(v_m).*n;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_K = ik(v,n,g_K,e_K)

% g_K = 36;
% e_K = -12.26;
I_K = g_K.*n.^4.*(v-e_K);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_Na = ina(v,m,h,g_Na,e_Na)

% g_Na = 120;
% e_Na = 117.56;
I_Na = g_Na.*m.^3.*h.*(v-e_Na);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_L = il(v, g_L, e_L)

% g_L = 0.3;
I_L = g_L.*(v-e_L);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%