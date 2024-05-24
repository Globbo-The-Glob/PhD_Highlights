%% Peano-Hasel Linear Series Elastic Actuator 
% Code models a single PH in series with a linear elastic element
% To do this, the code will redefine a new sub actuator and proceed to
% calculate the displacement for that actuator over a small zip step.
% At each step, it uses the previous displacement to calculate the elastic
% element tension force, asking if this is enough to stop deformation of
% the actuator.

%% TODO: 
% > [] Tidy up
% > [x] Apparently working
% > [x] Confirm sub_actuator method produces same results as OG_PH, get all
% plots
% > [x] Allow force to trail off
% > [] Automatic Truncation
% > [] Add Linear Elastic Element
% > [] Add scaling as one PH will not work with an
% > [] Model new behaviour
% > [] Expand to non-linear elastic element

%% Currently busted shit 

%% Housekeeping

clear 
close all 
clc

%% Preallocation for loop
Le = zeros(1,1000);
Lp = zeros(1,1000);
alpha0 = zeros(1,1000);
alpha = zeros(1,1000);
A = zeros(1,1000);
h = zeros(1,1000);

F = zeros(1,100000);
P = zeros(1,100000);
X = zeros(1,100000);

le = zeros(1,1000);
lp = zeros(1,1000);



%% Model Constants
Lp(1) = 0.02; % [m] Initial Pouch Length - From original kellaris paper (2018)
perm0 = 8.85418782e-12; %#[m-3kg-1s4A2]
permr = 2.2; % for BOPP
Eb = 700e6; %[V/m] the breakdown strength of BOPP
w = 0.01; %[m] from original kellaris paper
t = 18e-6; %[m] Bopp thickness used in initial designs
% k = 137e6; %[Nm-1] Achilles tendon average sourced from Litchwark 2005
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
Fb = (lambda)*(cos(alpha0(1))/(1-cos(alpha0(1)))); % Gives max actuator output in [N]

%% Sub-actuator calculations
% Assumes that each stage of actuator zipping can be defined as a new
% actuator with properties of unzippied region.
% Increase the le iteration to get higher resolution data
% --------------------------------------------------------------------
Comp = [0;0];
    
    subx = zeros(1,1000);
    x = zeros(1,100000);
    f = zeros(1,100000);
    p = zeros(1,100000);
    I = ones(1,20);
j=0
for k = 250:1000:10000
    j = j+1
    K(j) = k;
    
    for i  = 1:2
        I(i) = i;
        disp(i)
          
                end
            % Adds this loops vector to larger matrix

            F =[F; f];
            X =[X; x];
            P =[P; p];
        end
    end
end
F(1,:)=[];
F(:,1)=[];
P(1,:)=[];
P(:,1)=[];
X(1,:)=[];
X(:,1)=[];

Eng_Spring = 0.5.*K'.*X.^2;


strain = X(1,:)./h(1)*100; % Strain for actuator overall.
nF = F(1,:)./lambda;


figure

plot(strain,F(1,:))
hold on
for i = 2:size(K)
    disp(i)
    strain = X(i,:)./h(1)*100; % Strain for actuator overall.
    nF = F(i,:)./lambda;
    plot(strain,F(i,:))
end
hold off
title('Strain/Force SubActuator Model')
%xlim([0 25])
%ylim([0 15])
xlabel('Strain [%]')
ylabel('Force Output [N]')
