close all
clear

% Rigid rollrate solutions with symbolic computation

%% Data definition

run Data.m      % generic data input 
syms y_m        % definition of symbolic variables 

% Input dynamic pressure
fileID = fopen('div_pressure_Pa.txt','r');
if fileID == -1
    error('No input dynamic pressure found!');      % run Divergence.m before
else
    q_div_Pa = fscanf(fileID,'%f');                 % read dynamic divergence pressure in output file of 'Divergence.m'
    fclose(fileID);
end

% Evaluate roll dynamic pressure
q_roll_Pa = q_div_Pa * Roll_dynP_fraction;
Uinf_ms = sqrt((q_roll_Pa * 2) / rho_kgm3);


%% System matrix definition

%Ka
Ka_rigid =  - int ( cw_m .* CLalfaW .* (y_m).^2 .* (1 / Uinf_ms) , 0, bw_m) ...
            - int ( cc_m .* CLalfaC .* (y_m).^2 .* (1 / Uinf_ms) , 0, bc_m) ...
            - int ( cw_m .* CLalfaW .* (y_m).^2 .* (1 / Uinf_ms) , -bw_m, 0) ...
            - int ( cc_m .* CLalfaC .* (y_m).^2 .* (1 / Uinf_ms) , -bc_m, 0);

%f
f_p_rigid = + q_roll_Pa * int( cw_m .* CLbetaA .* Beta_rad .* y_m, ba_m, bw_m) ...
            + q_roll_Pa * int( cw_m .* CLbetaA .* (-Beta_rad) .* y_m, -bw_m, -ba_m) ;

H_rigid = -(q_roll_Pa*Ka_rigid);


%% Rollrate computation

res_rigid = H_rigid \ f_p_rigid;        % Linear system solution
p_rigid_rads = double(res_rigid);       % Rigid rollrate value [rad/s]
fprintf('Rigid rollrate is: %f rad/s\n',p_rigid_rads);