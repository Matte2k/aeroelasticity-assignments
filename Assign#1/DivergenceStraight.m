clear

% Divergence approximate Ritz solutions with symbolic computation for equivalent straight wing

%% Data definition

run Data_straight.m     % generic data input
syms y_mark_m           % definition of symbolic variables 
run ShapeFunctions.m    % shape function input


%% System matrix definition

%Ks
Ks_tt_w1 = + int(dN_theta' * GJw_Nm2 * dN_theta , 0, bw_mark_m);
Ks_tt_w2 = + int(dN_theta' * GJw_Nm2 * dN_theta , -bw_mark_m, 0);
Ks_zz_w1 = + int(d2N_z' * EIw_Nm2 * d2N_z , 0, bw_mark_m);
Ks_zz_w2 = + int(d2N_z' * EIw_Nm2 * d2N_z , -bw_mark_m, 0);

%Ka
Ka_tt_w1 = + int( (N_theta' * N_theta) * ew_m * cw_m * CLalfaW * (cos(SweepAng_rad_w1))^3 , 0,  bw_mark_m) ; 
Ka_tt_w2 = + int( (N_theta' * N_theta) * ew_m * cw_m * CLalfaW * (cos(SweepAng_rad_w2))^3 , -bw_mark_m, 0) ;

Ka_tz_w1 = - int( (N_theta' * dN_z) * ew_m * cw_m * CLalfaW * (cos(SweepAng_rad_w1)).^2 * sin(SweepAng_rad_w1), 0, bw_mark_m);
Ka_tz_w2 = - int( (N_theta' * dN_z) * ew_m * cw_m * CLalfaW * (cos(SweepAng_rad_w2)).^2 * sin(SweepAng_rad_w2), -bw_mark_m, 0);

Ka_zt_w1 =  + int( (N_z' * N_theta) * cw_m .* CLalfaW * (cos(SweepAng_rad_w1)).^2 , 0, bw_mark_m) ...
            - int( (dN_z' * N_theta) * ew_m * cw_m * CLalfaW .* (cos(SweepAng_rad_w1)).^2 .* sin(SweepAng_rad_w1) , 0, bw_mark_m);
Ka_zt_w2 =  + int( (N_z' * N_theta) * cw_m .* CLalfaW * (cos(SweepAng_rad_w2)).^2 , -bw_mark_m, 0) ...
            - int( (dN_z' * N_theta) * ew_m * cw_m * CLalfaW .* (cos(SweepAng_rad_w2)).^2 .* sin(SweepAng_rad_w2) , -bw_mark_m, 0);

Ka_zz_w1 =  - int( (N_z' * dN_z) * cw_m .* CLalfaW * sin(SweepAng_rad_w1)  .*(cos(SweepAng_rad_w1)) , 0, bw_mark_m) ...
            + int( (dN_z' * dN_z) * ew_m .* cw_m * CLalfaW * (sin(SweepAng_rad_w1))^2   .*(cos(SweepAng_rad_w1)) , 0, bw_mark_m);
Ka_zz_w2 =  - int( (N_z' * dN_z) * cw_m .* CLalfaW * sin(SweepAng_rad_w2)  .*(cos(SweepAng_rad_w2)) , -bw_mark_m, 0) ...
            + int( (dN_z' * dN_z) * ew_m .* cw_m * CLalfaW * (sin(SweepAng_rad_w2))^2   .*(cos(SweepAng_rad_w2)) , -bw_mark_m, 0);

% Matrix assembly
Ks = [  (Ks_tt_w1+Ks_tt_w2),        zeros(dimN_theta, dimN_z);  ...
        zeros(dimN_z,dimN_theta),   (Ks_zz_w1+Ks_zz_w2)      ];

Ka = [ (Ka_tt_w1+Ka_tt_w2),  (Ka_tz_w1+Ka_tz_w2); ...
       (Ka_zt_w1+Ka_zt_w2) , (Ka_zz_w1+Ka_zz_w2)];

Ka_disc=double(Ka);
Ks_disc=double(Ks);


%% Divergence pressure and velocity
[V, dynPd_Pa] = eig(Ks_disc, Ka_disc);                  % generalize eigenvalue problem for divergence
Ud_ms_vect = sqrt( (diag(dynPd_Pa)*2) / (rho_kgm3) );   % need to check for the minimun real eigenvalue


%% Output file

for i=1:(dimN_theta+dimN_z)             % elimination of complex eigenvalues
    if imag(Ud_ms_vect(i))~=0
            Ud_ms_vect(i)=inf;
    end
end

[Ud_ms, Ud_index] = min(Ud_ms_vect);    % first divergence speed value
Pd_Pa = dynPd_Pa(Ud_index,Ud_index);    % first divergence pressure value
fprintf('First divergence speed for equivalent straight wing is: %f m/s\n',Ud_ms);
fprintf('First divergence dynamic pressure for equivalent straight wing is: %f Pa\n',Pd_Pa);

% write straight divergence speed in output file
fileID = fopen('straight_div_speed_ms.txt','w');
text_output = '%f'; 
fprintf(fileID,text_output,Ud_ms);
fclose(fileID);

