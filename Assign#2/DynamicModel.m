close all
clear

% Dynamic system model

%% Data definition

run Data2.m             % generic data input
syms y_mark_m x_m       % definition of symbolic variables 
run ShapeFunctions_2.m    % shape function input

PTcstrC_m(1) = norm(PTcstrC_m(1));
CGw_m = norm(CGw_m);

%% Point evaluated shape
phi_wing = matlabFunction(N_phi);
phi_wing = phi_wing(0);
%phi_wing = N_phi;      % rigid fuselage

phi_canard = matlabFunction(N_phi);
phi_canard = phi_canard(PTcstrC_m(1));
%phi_canard = N_phi;    % rigid fuselage

theta_fuselage = matlabFunction(N_theta);
theta_fuselage = theta_fuselage(0);

z_fuselage = matlabFunction(N_z);
z_fuselage = z_fuselage(0);


%% Ks
Ks_tt_w1 = + int(dN_theta' * GJw_Nm2 * dN_theta , 0, bw_mark_m);
Ks_tt_w2 = + int(dN_theta' * GJw_Nm2 * dN_theta , -bw_mark_m, 0);
Ks_zz_w1 = + int(d2N_z' * EIw_Nm2 * d2N_z , 0, bw_mark_m);
Ks_zz_w2 = + int(d2N_z' * EIw_Nm2 * d2N_z , -bw_mark_m, 0);
Ks_pp = + int(dN_phi' * GJf_Nm2 * dN_phi , -PTcstrC_m(1,1), 0);
%Ks_pp = 0;     % rigid fuselage

% Matrix assembly
Ks = [  (Ks_tt_w1+Ks_tt_w2),        zeros(dimN_theta, dimN_z),      zeros(dimN_theta, dimN_phi); ...
        zeros(dimN_z,dimN_theta),   (Ks_zz_w1+Ks_zz_w2),            zeros(dimN_z,dimN_phi); ...
        zeros(dimN_phi,dimN_theta), zeros(dimN_phi, dimN_z),        Ks_pp ];

Ks_disc = double(Ks);


%% M

%M_tt
M_tt_w1 = double(   + int(N_theta' * N_theta * CGw_mark_m^2 * mw_kgm , 0, bw_mark_m) ...
                    - int(N_theta' * N_theta * CGw_mark_m * y_mark_m * mw_kgm * tan(SweepAng_rad_w1) , 0, bw_mark_m) * 2  ...      % 2 from math
                    + int(N_theta' * N_theta * y_mark_m^2 * mw_kgm * tan(SweepAng_rad_w1)^2 , 0, bw_mark_m) *1 ...     %check
                    + int(N_theta' * N_theta * Iw_kgm2, 0, bw_mark_m) );

M_tt_w2 = double(   + int(N_theta' * N_theta * CGw_mark_m^2 * mw_kgm , -bw_mark_m, 0) ...
                    - int(N_theta' * N_theta * CGw_mark_m * y_mark_m * mw_kgm * tan(SweepAng_rad_w2) , -bw_mark_m, 0) * 2  ...      % 2 from math
                    + int(N_theta' * N_theta * y_mark_m^2 * mw_kgm * tan(SweepAng_rad_w2)^2 , -bw_mark_m, 0) *1 ...      %check
                    + int(N_theta' * N_theta * Iw_kgm2, -bw_mark_m, 0) );

M_tt_f = double(   int(theta_fuselage' * theta_fuselage * (-x_m).^2 * mf_kgm * cos(SweepAng_rad)^2 , -lf_m, 0)  ...
                   + (theta_fuselage' * theta_fuselage * PTcstrC_m(1,1).^2 * Mc_kg * cos(SweepAng_rad)^2) );

M_tt = M_tt_w1 + M_tt_w2 + M_tt_f;


%M_tz (trasposto di zt)
M_tz_w1 = double(   - int(N_theta' * N_z * mw_kgm * CGw_mark_m, 0, bw_mark_m) ...
                    + int(N_theta' * N_z * mw_kgm * y_mark_m * tan(SweepAng_rad_w1), 0, bw_mark_m) *1  ...      %check
                    - int(N_theta' * dN_z * mw_kgm * y_mark_m * CGw_mark_m, 0, bw_mark_m) ...
                    + int(N_theta' * dN_z * mw_kgm * y_mark_m^2 * tan(SweepAng_rad_w1), 0, bw_mark_m) *1 );     %check

M_tz_w2 = double(   - int(N_theta' * N_z * mw_kgm * CGw_mark_m, 0, bw_mark_m) ...
                    + int(N_theta' * N_z * mw_kgm * y_mark_m * tan(SweepAng_rad_w2), -bw_mark_m, 0) *1 ...      %check
                    - int(N_theta' * dN_z * mw_kgm * y_mark_m * CGw_mark_m, -bw_mark_m, 0) ...
                    + int(N_theta' * dN_z * mw_kgm * y_mark_m^2 * tan(SweepAng_rad_w2), -bw_mark_m, 0) *1 );    %check

M_tz_f = double(   int(theta_fuselage' * z_fuselage * (-x_m) * mf_kgm * cos(SweepAng_rad), -lf_m, 0 ) ...
                   + (theta_fuselage' * z_fuselage * PTcstrC_m(1,1) * Mc_kg * cos(SweepAng_rad)) );

M_tz = M_tz_w1 + M_tz_w2 + M_tz_f;


%M_tp (trasposto di pt)
M_tp_w1 = double(   - int(N_theta' * phi_wing * mw_kgm * CGw_mark_m * y_mark_m * 1/cos(SweepAng_rad_w1) , 0, bw_mark_m) ...
                    + int(N_theta' * phi_wing * mw_kgm * y_mark_m^2 * 1/cos(SweepAng_rad_w1) * tan(SweepAng_rad_w1), 0, bw_mark_m) * 1 );  %check

M_tp_w2 = double(   - int(N_theta' * phi_wing * mw_kgm * CGw_mark_m * y_mark_m * 1/cos(SweepAng_rad_w2) , -bw_mark_m, 0) ...
                    + int(N_theta' * phi_wing * mw_kgm * y_mark_m^2 * 1/cos(SweepAng_rad_w2) * tan(SweepAng_rad_w2), -bw_mark_m, 0) * 1 );  %check

M_tp_f = double(   0 );

M_tp = M_tp_w1 + M_tp_w2 + M_tp_f;


%M_zz
M_zz_w1 = double(   + int(N_z' * N_z * mw_kgm , 0, bw_mark_m) ...
                    + int(N_z' * dN_z * mw_kgm * y_mark_m , 0, bw_mark_m)...
                    + int(dN_z' * N_z * mw_kgm * y_mark_m, 0, bw_mark_m) ...
                    + int(dN_z' * dN_z * mw_kgm * y_mark_m^2, 0, bw_mark_m) );

M_zz_w2 = double(   + int(N_z' * N_z * mw_kgm , -bw_mark_m, 0) ...
                    + int(N_z' * dN_z * mw_kgm * y_mark_m , -bw_mark_m, 0)...
                    + int(dN_z' * N_z * mw_kgm * y_mark_m, -bw_mark_m, 0) ...
                    + int(dN_z' * dN_z * mw_kgm * y_mark_m^2, -bw_mark_m, 0) );

M_zz_f = double(    (z_fuselage' * z_fuselage * mf_kgm * lf_m) + (z_fuselage' * z_fuselage * Mc_kg) );

M_zz = M_zz_w1 + M_zz_w2 + M_zz_f;


%M_zt  (trasposto di tz)
M_zt = M_tz';


%M_zp (trasposto di pz)
M_zp_w1 = double(   + int(N_z' * phi_wing * mw_kgm * y_mark_m * 1/cos(SweepAng_rad_w1), 0, bw_mark_m) ...
                    + int(dN_z' * phi_wing * mw_kgm * y_mark_m^2 * 1/cos(SweepAng_rad_w1), 0, bw_mark_m) );

M_zp_w2 = double(   + int(N_z' * phi_wing * mw_kgm * y_mark_m * cos(SweepAng_rad_w2), -bw_mark_m, 0) ...
                    + int(dN_z' * phi_wing * mw_kgm * y_mark_m^2 * 1/cos(SweepAng_rad_w2), -bw_mark_m, 0) );

M_zp_f = double(    0 );

M_zp = M_zp_w1 + M_zp_w2 + M_zp_f;


%M_pz (trasposto di zp)
M_pz = M_zp';


%M_pt (trasposto di tp)
M_pt = M_tp';


%M_pp
M_pp_w1 = double(   + int(phi_wing' * phi_wing * mw_kgm * y_mark_m.^2 * 1/cos(SweepAng_rad_w1)^2, 0, bw_mark_m) * 1 );  %check

M_pp_w2 = double(   + int(phi_wing' * phi_wing * mw_kgm * y_mark_m.^2 * 1/cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) * 1 );  %check

M_pp_f = double(    + int(N_phi' * N_phi * If_kgm2, -PTcstrC_m(1,1), 0)  );

M_pp = M_pp_w1 + M_pp_w2 + M_pp_f;


% Matrix assembly
M =   [  (M_tt),    (M_tz),    (M_tp); ...
         (M_zt),    (M_zz),    (M_zp); ...
         (M_pt),    (M_pz),    (M_pp)];

M_disc = double(M);


%% Ka
%Ka_z
Ka_zz_w1 = double(   - int(N_z' * dN_z * CLalfaW * cw_m * sin(SweepAng_rad_w1) * cos(SweepAng_rad_w1), 0, bw_mark_m) ...
                     + int(dN_z' * N_z * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w1)^2 * cos(SweepAng_rad_w1), 0, bw_mark_m) );

Ka_zz_w2 = double(   - int(N_z' * dN_z * CLalfaW * cw_m * sin(SweepAng_rad_w2) * cos(SweepAng_rad_w2), -bw_mark_m, 0) ...
                     + int(dN_z' * N_z * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w2)^2* cos(SweepAng_rad_w2), -bw_mark_m, 0) );

Ka_zz = Ka_zz_w1 + Ka_zz_w2;


Ka_zt_w1 = double(   + int(N_z' * N_theta * CLalfaW * cw_m * cos(SweepAng_rad_w1)^2, 0, bw_mark_m) ...
                     - int(dN_z' * N_theta * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w1) * cos(SweepAng_rad_w1)^2, 0, bw_mark_m) ...
                     + (z_fuselage' * theta_fuselage * CLalfaC * cc_m * cos(SweepAng_rad_w1) * bc_m) );

Ka_zt_w2 = double(   + int(N_z' * N_theta * CLalfaW * cw_m * cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) ...
                     - int(dN_z' * N_theta * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w2) * cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) ...
                     + (z_fuselage' * theta_fuselage * CLalfaC * cc_m * cos(SweepAng_rad_w2) * bc_m) );

Ka_zt = Ka_zt_w1 + Ka_zt_w2;


Ka_zp_w1 = 0;
Ka_zp_w2 = 0;
Ka_zp = (Ka_zp_w1 + Ka_zp_w2) * zeros(dimN_z,dimN_phi);


%Ka_t
Ka_tz_w1 = double(   - int(N_theta' * dN_z * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w1) * cos(SweepAng_rad_w1)^2, 0, bw_mark_m) );

Ka_tz_w2 = double(   - int(N_theta' * dN_z * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w2) * cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) );

Ka_tz = Ka_tz_w1 + Ka_tz_w2;


Ka_tt_w1 = double(   + int(N_theta' * N_theta * CLalfaW * cw_m * ew_m * cos(SweepAng_rad_w1)^3, 0, bw_mark_m ) ...
                     + (theta_fuselage' * theta_fuselage * CLalfaC * cc_m * PTcstrC_m(1,1) * cos(SweepAng_rad_w1)^2 * bc_m) );

Ka_tt_w2 = double(   + int(N_theta' * N_theta * CLalfaW * cw_m * ew_m * cos(SweepAng_rad_w2)^3, -bw_mark_m, 0) ...
                     + (theta_fuselage' * theta_fuselage * CLalfaC * cc_m * PTcstrC_m(1,1) * cos(SweepAng_rad_w1)^2 * bc_m) );

Ka_tt = Ka_tt_w1 + Ka_tt_w2;


Ka_tp_w1 = 0;
Ka_tp_w2 = 0;
Ka_tp = (Ka_tp_w1 + Ka_tp_w2) * zeros(dimN_theta,dimN_phi);


%Ka_p
Ka_pz_w1 = double(   - int(phi_wing' * dN_z * CLalfaW * cw_m * y_mark_m * sin(SweepAng_rad_w1) * cos(SweepAng_rad_w1)^2, 0, bw_mark_m) );
Ka_pz_w2 = double(   - int(phi_wing' * dN_z * CLalfaW * cw_m * y_mark_m * sin(SweepAng_rad_w2) * cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) );
Ka_pz = Ka_pz_w1 + Ka_pz_w2;


Ka_pt_w1 = double(   + int(phi_wing' * N_theta * CLalfaW * cw_m * y_mark_m * cos(SweepAng_rad_w1)^3, 0, bw_mark_m) ...
                     + int(phi_canard' * theta_fuselage * cc_m * CLalfaC * y_mark_m * cos(SweepAng_rad_w1), 0, bw_mark_m) );
            
Ka_pt_w2 = double(   + int(phi_wing' * N_theta * CLalfaW * cw_m * y_mark_m * cos(SweepAng_rad_w2)^3, -bw_mark_m, 0) ...
                     + int(phi_canard' * theta_fuselage * cc_m * CLalfaC * y_mark_m * cos(SweepAng_rad_w2), -bw_mark_m, 0) );

Ka_pt = Ka_pt_w1 + Ka_pt_w2;


Ka_pp_w1 = 0;
Ka_pp_w2 = 0;
Ka_pp = (Ka_pp_w1 + Ka_pp_w2) * zeros(dimN_phi,dimN_phi);

% Matrix assembly
Ka =   [ Ka_tt,    Ka_tz,  Ka_tp; ...
         Ka_zt,    Ka_zz,  Ka_zp; ...
         Ka_pt,    Ka_pz,  Ka_pp ];

Ka_disc = double(Ka);


%% Ca matrix definition

%Ca_z
Ca_zz_w1 = double(  - int(N_z' * N_z * CLalfaW * cw_m * cos(SweepAng_rad_w1), 0, bw_mark_m) ...
                    + int(dN_z' * N_z * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w1) * cos(SweepAng_rad_w1), 0, bw_mark_m) ...
                    - (z_fuselage' * z_fuselage * CLalfaC * cc_m * bc_m) );
Ca_zz_w2 = double(  - int(N_z' * N_z * CLalfaW * cw_m * cos(SweepAng_rad_w2), -bw_mark_m, 0) ...
                    + int(dN_z' * N_z * CLalfaW * cw_m * ew_m * sin(SweepAng_rad_w2) * cos(SweepAng_rad_w2), -bw_mark_m, 0) ...
                    - (z_fuselage' * z_fuselage * CLalfaC * cc_m * bc_m) );
Ca_zz = Ca_zz_w1 + Ca_zz_w2;

Ca_zt = zeros(dimN_z,dimN_theta);
Ca_zp = zeros(dimN_z,dimN_phi);

%Ca_t
Ca_tz_w1 = double(  - int(N_theta' * N_z * CLalfaW * cw_m * ew_m * cos(SweepAng_rad_w1)^2, 0, bw_mark_m) ...
                    - (theta_fuselage' * z_fuselage * CLalfaC * cc_m * PTcstrC_m(1,1) * cos(SweepAng_rad_w1) * bc_m ) );
Ca_tz_w2 = double(  - int(N_z' * N_theta * CLalfaW * cw_m * ew_m * cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) ...
                    - (theta_fuselage' * z_fuselage * CLalfaC * cc_m * PTcstrC_m(1,1) * cos(SweepAng_rad_w2) * bc_m ) );
Ca_tz = Ca_tz_w1 + Ca_tz_w2;

Ca_tt =  zeros(dimN_theta,dimN_theta);
Ca_tp = zeros(dimN_theta,dimN_phi);

%Ca_p
Ca_pz_w1 = double(  - int(phi_wing' * N_z* CLalfaW * cw_m * y_mark_m * cos(SweepAng_rad_w1)^2, 0, bw_mark_m) ...
                    - int(phi_canard' * z_fuselage * CLalfaC * cc_m * y_mark_m, 0, bc_m) );
Ca_pz_w2 = double(  - int(phi_wing' * N_z * CLalfaW * cw_m * y_mark_m * cos(SweepAng_rad_w2)^2, -bw_mark_m, 0) ...
                    - int(phi_canard' * z_fuselage * CLalfaC * cc_m * y_mark_m, -bc_m, 0) );
Ca_pz = Ca_pz_w1 + Ca_pz_w2;

Ca_pt = zeros(dimN_phi,dimN_theta);    
Ca_pp = zeros(dimN_phi,dimN_phi);


% Matrix assemble
Ca =   [ Ca_tt,    Ca_tz,  Ca_tp; ...
         Ca_zt,    Ca_zz,  Ca_zp; ...
         Ca_pt,    Ca_pz,  Ca_pp ];

Ca_disc = double(Ca);

