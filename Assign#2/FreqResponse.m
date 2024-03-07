% Frequency response computation for system:  M*q'' + C*q' + K*q = f(w)

%% Speed and frequency conditions

U_resp_ms = 30;                     % speed at which compute frequency response [m/s]
q_resp_Pa = (1/2 * rho_kgm3 * U_resp_ms^2);

InputFreq_Hz = x1 * 15;             % input frequency [Hz]
alpha_canard_rad = deg2rad(5);      % canard deflection [rad]

%% Forcing term assembly

%F
F_theta =   double(     theta_fuselage' * cos(SweepAng_rad_w1) * PTcstrC_m(1,1) * cc_m * CLalfaC * alpha_canard_rad * q_resp_Pa * bc_m) ...
            + double(   theta_fuselage' * cos(SweepAng_rad_w2) * PTcstrC_m(1,1) * cc_m * CLalfaC * alpha_canard_rad * q_resp_Pa * bc_m);

F_z =   double(     z_fuselage' * cc_m * CLalfaC * alpha_canard_rad * q_resp_Pa * bc_m) ...
        + double(   z_fuselage' * cc_m * CLalfaC * alpha_canard_rad * q_resp_Pa * bc_m);

F_phi =   double(   int(phi_canard' * y_mark_m * cc_m * CLalfaC * alpha_canard_rad * q_resp_Pa, 0, bw_mark_m)) ...
          + double( int(phi_canard' * y_mark_m * cc_m * CLalfaC * alpha_canard_rad * q_resp_Pa, -bw_mark_m, 0));

F_disc = [F_theta; F_z; F_phi];


%% Modal semplification of the system

Minv_modal = inv(M_modal);      % inverse of modal mass
I_modal = eye(Sys_dim);         % identity matrix for state space system
SS_Sys_dim = 2 * Sys_dim;       % modal system element dimension

Cs_modal = diag(2 .* diag(M_modal) .* damp_fraction .* EigFreqReal);        % modal damping matrix  
Ca_modal = EigVect_ort' * Ca_disc * EigVect_ort;                            % modal aerodynamic damping matrix
C_modal = (Cs_modal) - ( 1/2 * rho_kgm3 * U_resp_ms) * (Ca_modal)  ;        % modal total damping

K_modal = EigVect_ort' * (Ks_disc - (1/2 * rho_kgm3 * U_resp_ms^2) * (Ka_disc)) * EigVect_ort;    % modal total stiffness

F_modal = EigVect_ort' * F_disc;    % modal forcing term

T = diag(F_modal);

%% Modal state space system definition

A = [   zeros(Sys_dim),         I_modal; ...
        (-Minv_modal)*K_modal,  (-Minv_modal)*C_modal ];

B = [   zeros(Sys_dim);         (Minv_modal)*T  ];

C = [   ones(1,Sys_dim),        zeros(1,Sys_dim)      ];    % looking for state

D = [   zeros(1,Sys_dim)    ];

StateSpaceSys_time = ss(A,B,C,D);                           % dynamic state space sys
StateSpaceSys_frq  = frd(StateSpaceSys_time,InputFreq_Hz);  % frequency-response sys
[Sys_mag,Sys_phase,Sys_wout] = bode(StateSpaceSys_frq);     % magnitude and phase computation

% bodes plots
freq_vect_Hz = logspace(-2,2,50);
StateSpaceSys_frq_disc  = frd(StateSpaceSys_time,freq_vect_Hz);
bode(StateSpaceSys_frq_disc)


%% Text output

% result print in command window
fprintf('System pitch (z_theta) frequencies response at %f Hz are:\n', Sys_wout);
fprintf('Magnitudes :  ');
fprintf('%f   ',Sys_mag(1:dimN_theta));
fprintf('\nPhases     :  ');
fprintf('%f   ',Sys_phase(1:dimN_theta));
fprintf('\n\n');

fprintf('System plunge (z_z) frequencies response at %f Hz are:\n', Sys_wout);
fprintf('Magnitudes :  ');
fprintf('%f   ',Sys_mag((dimN_theta+1):(dimN_theta+dimN_z)));
fprintf('\nPhases     :  ');
fprintf('%f   ',Sys_phase((dimN_theta+1):(dimN_theta+dimN_z)));
fprintf('\n\n');

fprintf('System roll (z_phi) frequencies response at %f Hz are:\n', Sys_wout);
fprintf('Magnitudes :  ');
fprintf('%f   ',Sys_mag((dimN_theta+dimN_z+1):(dimN_theta+dimN_z+dimN_phi)));
fprintf('\nPhases     :  ');
fprintf('%f   ',Sys_phase((dimN_theta+dimN_z+1):(dimN_theta+dimN_z+dimN_phi)));
fprintf('\n\n');
