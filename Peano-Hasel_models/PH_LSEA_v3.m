%% Peano-Hasel Linear Series Elastic Actuator 
% Code models a single PH in series with a linear elastic element
% To do this, the code will redefine a new sub actuator and proceed to
% calculate the displacement for that actuator over a small zip step.
% At each step, it uses the previous displacement to calculate the elastic
% element tension force, asking if this is enough to stop deformation of
% the actuator.
% TODO: 
% > Debug
% > Confirm sub_actuator method produces same results as OG_PH, get all
% plots
% > Add scaling as one PH will not work with spring
% > Model new behaviour

%% Currested busted shit 
% The eqn solver at the end of the while loop is fucked. Cause unknown.

%% Housekeeping

clear 
close all 
clc

%% Preallocation for loop
Le = ones(1,1000);
Lp = ones(1,1000);
alpha0 = ones(1,1000);
alpha = ones(1,1000);
A = ones(1,1000);
h = ones(1,1000);
F = ones(1,1000);
P = ones(1,1000);
le = ones(1,1000);
lp = ones(1,1000);
subx = zeros(1,1000);
x = zeros(1,1000);
%% Model Constants
Lp(1) = 0.02; % [m] Initial Pouch Length - From original kellaris paper (2018)
perm0 = 8.85418782e-12; %#[m-3kg-1s4A2]
permr = 2.2; % for BOPP
Eb = 700e6; %[V/m] the breakdown strength of BOPP
w = 0.01; %[m] from original kellaris paper
t = 18e-6; %[m] Bopp thickness used in initial designs
k = 137e6; %[Nm-1] Achilles tendon average sourced from Litchwark 2005
rho_f = 903; %[kgm-3] Bopp Density TEKRA (2020)
rho_o = 920; %[kgm-3] FR3 Dielectric oil Density Cargill (2018)
alphaf = pi/2; %[rad] Assume circular cross section on stroke completion
%% Swept constants
Eper = 0.5; %1x10 percentage of electrode coverage
% E = linspace(0,Eb,10); %1x10 Varying field strength [V/m]
E = 700e6;

%Combine for force equation 
lambda = E^2*perm0*permr*w*t;

%% Setup
Eper = Eper';
%Electrode length initial [m] 1x10
Le(1) = Lp(1).*Eper;

%% Solves fat optimal fill eqn for alpha0
alpeqn = ((1-Eper).^2)*(2/pi);
syms y
eqn = (y-sin(y)*cos(y))/y^2 == alpeqn; %Sets up alpha0 = y only eqn
alpha0(1) = double(vpasolve(eqn,y, [0,1.6])); % solve eqn for alpha0 in realistic interval. Double(): sym -> DubPre [rads] 10x1 or 1x10


%Height
h(1) = Lp(1)*(sin(alpha0(1))/alpha0(1)); % Height [m]
%Area w/ restraints 
A(1) = (Lp(1) - Le(1))^2 /pi; % Area [m2]

% Deformation setup
lp(1) = Lp(1); % no deformation in first stage
le(1) = 0; % no intial zipping
alpha(1) = alpha0(1);
%Voltage
V = 2*t*E ;


%% Max Force of actuator
Fb = (w*t*E^2*perm0*permr)*(cos(alpha0(1))/(1-cos(alpha0(1)))); % Gives max actuator output in [N]

%% Sub-actuator calculations
% Assumes that each stage of actuator zipping can be defined as a new
% actuator with properties of unzippied region.
% Increase the le iteration to get higher resolution data
% 
n = 1;
%while n < 1000
while n<10000
    n = n+1; % Iterates n
    
    Le(n) = Le(1) - le(n-1); % New electrode coverage is intial minus previous zip
    
    Lp(n) = lp(n-1); % New pouch length is previous final pouch length
    
    Eper(n) = Le(n)/Lp(n); % Calculates electrode coverage percentage
    
    alpha0(n) = alpha(n-1);
    
    A(n) = 0.5*Lp(n)^2*((alpha0(n) - sin(alpha0(n))*cos(alpha0(n)))/alpha0(n)^2); % New cross section area from non-optimum fill equation
    
    h(n) = Lp(n)*(sin(alpha0(n))/alpha0(n));
    
    F(n) = lambda*(cos(alpha0(n))/(1-cos(alpha0(n)))); % Calculates force generated by sub-actuator
    
    P(n) = 0.1*9.81; %x(n-1)*k; % Calculates current value of spring tension\ currently a pure weight
    
    if F(n) <= P(n) % If the actuator and spring are in eq, no displacement, loop ends
       
        break
        
    elseif Le == 0 % If the electrode coverage is now 0, there is no more zipping
   
        break
        
    else % If actuator is still out forcing spring, there is displacement,
        
        le(n) = le(n-1) + 0.00001*(Le(1)); % Incriments zip length by adding a little bit of the original electrode length
        
        lp(n) = Lp(n) - le(n); % New pouch length after zip iteration
        
        syms y                                                 
        lpeqn = sqrt((2*A(n)*y^2)/(y-sin(y)*cos(y))) == lp(n); % |Solves lp(alpha) eqn for new alpha 
        yans = double(vpasolve(lpeqn,y,[0,1.6]));

        alpha(n) = yans;
        
        subx(n) = h(n)-(lp(n)*(sin(alpha(n))/alpha(n)) + le(n)); % Displacement of sub-actuator 
        x(n) = x(n-1)+subx(n); % Total displacement of actuator.
    end
end    

strain = x./h(1)*100; % Strain for actuator overall.

figure
plot(strain,F)
