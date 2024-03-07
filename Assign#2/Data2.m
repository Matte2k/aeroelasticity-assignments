close all
clear

%% PERSONAL data
x1 = 1.2;             %constant one from personal code
x2 = 1.1;             %constant two from personal code

rho_kgm3 = 1.225;               % divercence speed density [kg/m3]
damp_fraction = 0.03;           % percentage of structural proportional damping

%% WING data
SweepAng_deg = -30;                             % wing sweep angle [deg]
SweepAng_rad = deg2rad(SweepAng_deg);           % wing sweep angle [rad]
SweepAng_rad_w1 = deg2rad(SweepAng_deg);        % wing1 sweep angle [rad]
SweepAng_rad_w2 = deg2rad(-SweepAng_deg);       % wing2 sweep angle [rad]

bw_m = 6.10;                            % wing semispan [m]
bw_mark_m = bw_m / cos(SweepAng_rad);

cw_m = 3.05;                            % wing chord [m]
cw_mark_m = cw_m * cos(SweepAng_rad);

PTcstrW_m = [0,0];                          % wing point of constraint (line vector x,y) [m]
CGw_mark_m = 0.42 * x2/cw_mark_m;           % wing CG position structural system [m]
CGw_m = CGw_mark_m / cos(SweepAng_rad);     % wing CG position [m]
mw_kgm = 400;                               % wing mass per unit span [kg/m]

EAw_m = 1/2 * cw_m;                     % wing elastic axis [m]
EIw_Nm2 = x1 * 4.5 * 1e6;               % wing bending stiffness [Nm2]
GJw_Nm2 = x2 * 7 * 1e6;                 % wing torsional stiffness [Nm2]
Iw_kgm2 = 320;                          % wing moment of inertia twist per unit span [kg m2]
Jw_kgm2 = mw_kgm * bw_mark_m^2 / 3;     % wing moment of inertia roll per unit span [kg m2]
CLalfaW = 2*pi;                         % wing CL_alfa

ACw_m = 1/4 * cw_m;                 % wing aerodinamic centre positin [m]
ew_m = EAw_m - ACw_m;               % wing distance between ElasticAxis and Aerodynamic centre [m]


%% FUSELAGE data
lf_m = 9.15;                       % fuselage length [m]
mf_kgm = x1 * 500;                 % fuselage mass per unit span [kg/m]

GJf_Nm2 = 12e6;                    % fuselage torsional stiffness [Nm2]
If_kgm2 = 3000;                    % fuselage moment of inertia per unit span [kg m2]


%% CANARD data
bc_m = 1.525;                       % canard semispan [m]
cc_m = 3.05;                        % canard chord [m]
PTcstrC_m = [-4.525,0];             % canard point of constraint (line vector x,y) [m]
CLalfaC = 2*pi;                     % canard CL_alfa
ACc_m = 1/4 * cc_m;                 % canard aerodinamic centre positin [m]
Mc_kg = 300;                        % canard total lumped mass applied in PTcstrC_m [kg]
