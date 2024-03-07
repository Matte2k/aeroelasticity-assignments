close all
clear 

% Rollrate approximate Ritz solutions with symbolic computation

%% Data definition
run Data.m              % generic data input
syms y_mark_m           % definition of symbolic variables
run ShapeFunctions.m    % shape function input

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
Uinf_ms = sqrt((q_roll_Pa * 2)/1.225);
Uinf_mark_ms = Uinf_ms * cos(SweepAng_rad);
 

%% System matrix definition

%Ks
Ks_tt_w1 =  int(GJw_Nm2 .* (dN_theta' * dN_theta) , y_mark_m, 0, bw_mark_m);
Ks_tt_w2 =  int(GJw_Nm2 .* (dN_theta' * dN_theta) , y_mark_m, -bw_mark_m, 0) ;

Ks_zz_w1 =  int(EIw_Nm2 .* (d2N_z' * d2N_z) , y_mark_m, 0, bw_mark_m);
Ks_zz_w2 =  int(EIw_Nm2 .* (d2N_z' * d2N_z) , y_mark_m, -bw_mark_m, 0);


%Ka
Ka_tt_w1 = + int( (N_theta' * N_theta) .* ew_m .* cw_m .* CLalfaW .* (cos(SweepAng_rad_w1)).^3 , y_mark_m, 0,  bw_mark_m);                                          
Ka_tt_w2 = + int( (N_theta' * N_theta) .* ew_m .* cw_m .* CLalfaW .* (cos(SweepAng_rad_w2)).^3 , y_mark_m, -bw_mark_m,  0);

Ka_tz_w1 = - int( (N_theta' * dN_z) .* ew_m .* cw_m * CLalfaW * (cos(SweepAng_rad_w1)).^2 .* sin(SweepAng_rad_w1) , y_mark_m, 0, bw_mark_m);                        
Ka_tz_w2 = - int( (N_theta' * dN_z) .* ew_m .* cw_m * CLalfaW * (cos(SweepAng_rad_w2)).^2 .* sin(SweepAng_rad_w2) , y_mark_m, -bw_mark_m, 0);                                 

Ka_tp_w1 =  - int ( N_theta' .* cw_m .* ew_m .* CLalfaW .* y_mark_m .* (1 / Uinf_ms) .* cos(SweepAng_rad_w1).^3, y_mark_m, 0, bw_mark_m);                           
Ka_tp_w2 =  - int ( N_theta' .* cw_m .* ew_m .* CLalfaW .* y_mark_m .* (1 / Uinf_ms) .* cos(SweepAng_rad_w2).^3, y_mark_m, -bw_mark_m, 0);                       

Ka_zt_w1 =  + int( (N_z' * N_theta) .* cw_m .* CLalfaW .* (cos(SweepAng_rad_w1)).^2 , y_mark_m, 0, bw_mark_m) ...                                                   
            - int( (dN_z' * N_theta) .* ew_m .* cw_m .* CLalfaW .* (cos(SweepAng_rad_w1)).^2 .* sin(SweepAng_rad_w1) , y_mark_m, 0, bw_mark_m);                     
Ka_zt_w2 =  + int( (N_z' * N_theta) .* cw_m .* CLalfaW .* (cos(SweepAng_rad_w2)).^2 , y_mark_m, -bw_mark_m, 0) ...
            - int( (dN_z' * N_theta) .* ew_m .* cw_m .* CLalfaW .* (cos(SweepAng_rad_w2)).^2 .* sin(SweepAng_rad_w2) , y_mark_m, -bw_mark_m, 0);  

Ka_zz_w1 =  - int( (N_z' * dN_z) .* cw_m .* CLalfaW .* cos(SweepAng_rad_w1) .* sin(SweepAng_rad_w1) , y_mark_m, 0, bw_mark_m) ...                                   
            + int( (dN_z' * dN_z) .* ew_m .* cw_m .* CLalfaW .* cos(SweepAng_rad_w1) .* (sin(SweepAng_rad_w1)).^2 , y_mark_m, 0, bw_mark_m);                        
Ka_zz_w2 =  - int( (N_z' * dN_z) .* cw_m .* CLalfaW .* cos(SweepAng_rad_w2) .* sin(SweepAng_rad_w2) , y_mark_m, -bw_mark_m, 0) ... 
            + int( (dN_z' * dN_z) .* ew_m .* cw_m .* CLalfaW .* cos(SweepAng_rad_w2) .* (sin(SweepAng_rad_w2)).^2 , y_mark_m, -bw_mark_m, 0);  

Ka_zp_w1 =  - int ( N_z' .* cw_m .* CLalfaW .* y_mark_m .* (1 / Uinf_ms) .* (cos(SweepAng_rad_w1)).^2, y_mark_m, 0, bw_mark_m) ...                                  
            + int( (dN_z') .* ew_m .* cw_m .* CLalfaW .* y_mark_m .* (1 / Uinf_ms) .* (cos(SweepAng_rad_w1)).^2 .* sin(SweepAng_rad_w1), y_mark_m, 0, bw_mark_m);   
Ka_zp_w2 =  - int ( N_z' .* cw_m .* CLalfaW .* y_mark_m .* (1 / Uinf_ms) .* (cos(SweepAng_rad_w2)).^2, y_mark_m, -bw_mark_m, 0) ...                                     
            + int( (dN_z') .* ew_m .* cw_m .* CLalfaW .* y_mark_m .* (1 / Uinf_ms) .* (cos(SweepAng_rad_w2)).^2 .* sin(SweepAng_rad_w2), y_mark_m, -bw_mark_m, 0);      
    
Ka_pt_w1 =  + int ( N_theta .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w1).^3, y_mark_m, 0, bw_mark_m);
Ka_pt_w2 =  + int ( N_theta .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w2).^3, y_mark_m, -bw_mark_m, 0) ;                                                            

Ka_pz_w1 =  - int ( dN_z .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w1).^2 .* sin(SweepAng_rad_w1) , y_mark_m, 0, bw_mark_m);
Ka_pz_w2 =  - int ( dN_z .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w2).^2 .* sin(SweepAng_rad_w2) , y_mark_m, -bw_mark_m, 0);                                         

Ka_pp_w1 =  - int ( cw_m .* CLalfaW .* (y_mark_m).^2 .* (1 / Uinf_ms) .* cos(SweepAng_rad_w1).^3 , y_mark_m, 0, bw_mark_m) ...
            - int ( cc_m .* CLalfaC .* (y_mark_m).^2 .* (1 / Uinf_ms), y_mark_m, 0, bc_m);           
Ka_pp_w2 =  - int ( cw_m .* CLalfaW .* (y_mark_m).^2 .* (1 / Uinf_ms) .* cos(SweepAng_rad_w2).^3 , y_mark_m, -bw_mark_m, 0) ...                                    
            - int ( cc_m .* CLalfaC .* (y_mark_m).^2 .* (1 / Uinf_ms), y_mark_m, -bc_m, 0);                                                                         


%f
f_t_w1 =    + int( (N_theta') .* cw_m .* ew_m .* CLbetaA .* (Beta_rad) .* cos(SweepAng_rad_w1).^2, y_mark_m, ba_mark_m, bw_mark_m) ...                          
            + int( (N_theta') .* cw_m.^2 .* CMbetaA .* (Beta_rad) .* cos(SweepAng_rad_w1).^2, y_mark_m, ba_mark_m, bw_mark_m);                                  
f_t_w2 =    + int( (N_theta') .* cw_m .* ew_m .* CLbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2).^2, y_mark_m, -bw_mark_m, -ba_mark_m) ...
            + int( (N_theta') .* cw_m.^2 .* CMbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2).^2, y_mark_m, -bw_mark_m, -ba_mark_m);

f_z_w1 =    + int( (N_z') .* cw_m .* CLbetaA .* (Beta_rad) .* cos(SweepAng_rad_w1), y_mark_m, ba_mark_m, bw_mark_m)...                                        
            - int( (dN_z') .* cw_m.^2 .* CMbetaA .* (Beta_rad) .* cos(SweepAng_rad_w1) .* sin(SweepAng_rad_w1), y_mark_m, ba_mark_m, bw_mark_m) ...
            - int( (dN_z') .* ew_m .* cw_m .* CLbetaA .* (Beta_rad) .* cos(SweepAng_rad_w1) .* sin(SweepAng_rad_w1), y_mark_m, ba_mark_m, bw_mark_m);
f_z_w2 =    + int( (N_z') .* cw_m .* CLbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2), y_mark_m, -bw_mark_m, -ba_mark_m)...
            - int( (dN_z') .* cw_m.^2 .* CMbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2) .* sin(SweepAng_rad_w2), y_mark_m, -bw_mark_m, -ba_mark_m) ...
            - int( (dN_z') .* ew_m .* cw_m .* CLbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2) .* sin(SweepAng_rad_w2), y_mark_m, -bw_mark_m, -ba_mark_m);

f_p_w1 = + int( cw_m .* CLbetaA .* (Beta_rad) .* cos(SweepAng_rad_w1).^2 .* y_mark_m, y_mark_m, ba_mark_m, bw_mark_m);
f_p_w2 = + int( cw_m .* CLbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2).^2 .* y_mark_m, y_mark_m, -bw_mark_m, -ba_mark_m);


% Matrix assembly
Ks = [(Ks_tt_w1 + Ks_tt_w2),        zeros(dimN_theta, dimN_z),      zeros(dimN_theta, 1); ...
      zeros(dimN_z,dimN_theta),     (Ks_zz_w1 + Ks_zz_w2),          zeros(dimN_z,1); ...
      zeros(1,dimN_theta),          zeros(1,dimN_z),                0];
  
Ka = [(Ka_tt_w1 + Ka_tt_w2),    (Ka_tz_w1 + Ka_tz_w2),  (Ka_tp_w1 + Ka_tp_w2); ...
      (Ka_zt_w1 + Ka_zt_w2),    (Ka_zz_w1 + Ka_zz_w2),  (Ka_zp_w1 + Ka_zp_w2); ...
      (Ka_pt_w1 + Ka_pt_w2),    (Ka_pz_w1 + Ka_pz_w2),  (Ka_pp_w1 + Ka_pp_w2)];

f = [(f_t_w1 + f_t_w2); (f_z_w1 + f_z_w2); (f_p_w1 + f_p_w2)] .* q_roll_Pa;
H = Ks - (q_roll_Pa * Ka);

% Matrix discretization
Ka_disc = double(Ka);
Ks_disc = double(Ks);
H_disc = double(H);
f_disc = double(f);


%% Rollrate computation
res = H_disc \ f_disc;                    % Linear system solution
p_rads = res((1+dimN_z+dimN_theta),1);    % Elastic rollrate value [rad/s] 
fprintf('Elastic rollrate is: %f rad/s\n\n',p_rads);
