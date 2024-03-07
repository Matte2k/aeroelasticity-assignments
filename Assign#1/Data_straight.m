clear
close all

%% PERSONAL data
x1 = 1.2;             %constant one from personal code
x2 = 1.1;             %constant two from personal code

rho_kgm3 = 1.225;               % divercence speed density
Roll_dynP_fraction = 0.85;      % percentage of dynPressure in roll

%% WING data
SweepAng_deg = 0;                               % wing sweep angle [deg]
SweepAng_rad = deg2rad(SweepAng_deg);           % wing sweep angle [rad]
SweepAng_rad_w1 = deg2rad(SweepAng_deg);        % wing1 sweep angle [rad]
SweepAng_rad_w2 = deg2rad(-SweepAng_deg);       % wing2 sweep angle [rad]

bw_m = 6.10;                        % wing semispan [m]
cw_m = 3.05;                        % wing chord [m]
PTcstrW_m = [0,0];                  % wing point of constraint (line vector x,y) [m]
EAw_m = 1/2 * cw_m;                 % wing elastic axis [m]
EIw_Nm2 = x1 * 4.5 * 1e6;           % wing bending stiffness [Nm2]
GJw_Nm2 = x2 * 7 * 1e6;             % wing torsional stiffness [Nm2]
CLalfaW = 2*pi;                     % wing CL_alfa
ACw_m = 1/4 * cw_m;                 % wing aerodinamic centre positin [m]
ew_m = EAw_m - ACw_m;               % wing distance between ElasticAxis and Aerodynamic centre [m]


% Projection on structural system
bw_mark_m = bw_m / cos(SweepAng_rad);
cw_mark_m = cw_m * cos(SweepAng_rad);
CLalfaW_mark = CLalfaW / cos(SweepAng_rad);
EAw_mark_m = 1/2 * cw_mark_m;
ACw_mark_m = 1/4 * cw_mark_m;
ew_mark_m = EAw_mark_m - ACw_mark_m;

%% FUSELAGE data
lf_m = 9.15;                       % fuselage length [m]
GJf_Nm2 = 12e6;                    % fuselage torsional stiffness [Nm2]

%% ALEIRON data
ba_m = 1/2 * bw_m;                                  % aleiron semispan [m]
ca_m = 1/4 * cw_m;                                  % aleiron chord [m]
Ea = ca_m/cw_m;                                     % constant E for aeileron
CLbetaA = 2*(acos(1-(2*Ea))+2*sqrt(Ea*(1-Ea)));     % aleiron CL_beta
CMbetaA = -2*(1-Ea)*sqrt(Ea*(1-Ea));                % aleiron Cm_beta
Beta_rad = deg2rad(1);                              % aleiron deflection [rad]

% Projection on structural system
ba_mark_m = 1/2 * bw_mark_m;
ca_mark_m = 1/4 * cw_mark_m;
Ea_mark = ca_mark_m / cw_mark_m;
CLbetaA_mark = 2*(acos(1-(2*Ea_mark))+2*sqrt(Ea_mark*(1-Ea_mark)));
CMbetaA_mark = -2*(1-Ea_mark)*sqrt(Ea_mark*(1-Ea_mark));

%% CANARD data
bc_m = 1.525;                       % canard semispan [m]
cc_m = 3.05;                        % canard chord [m]
PTcstrC_m = [-4.525,0];             % canard point of constraint (line vector x,y) [m]
CLalfaC = 2*pi;                     % canard CL_alfa
ACc_m = 1/4 * cc_m;                 % canard aerodinamic centre positin [m]   
