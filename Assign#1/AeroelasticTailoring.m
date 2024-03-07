clear all
% Approximate Ritz solutions with symbolic computation

%% Data definition
run Data.m              % generic data input
syms y_mark_m  Cw       % definition of symbolic variables
run ShapeFunctions.m    % shape function input

%% System matrix definition

%Ks
Ks_tt_w1 = + int(dN_theta' * GJw_Nm2 * dN_theta , 0, bw_mark_m);
Ks_tt_w2 = + int(dN_theta' * GJw_Nm2 * dN_theta , -bw_mark_m, 0);

Ks_tz_w1 = + int(dN_theta' * Cw * d2N_z , 0, bw_mark_m);
Ks_tz_w2 = + int(dN_theta' * Cw * d2N_z , -bw_mark_m, 0);

Ks_zt_w1 = + int(d2N_z' * Cw * dN_theta , 0, bw_mark_m);
Ks_zt_w2 = + int(d2N_z' * Cw * dN_theta , -bw_mark_m, 0);

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

%Build matrix
Ks = [ (Ks_tt_w1+Ks_tt_w2),   (Ks_tz_w1+Ks_tz_w2);  ...
       (Ks_zt_w1+Ks_zt_w2),   (Ks_zz_w1+Ks_zz_w2)];

Ka = [ (Ka_tt_w1+Ka_tt_w2),  (Ka_tz_w1+Ka_tz_w2); ...
       (Ka_zt_w1+Ka_zt_w2) , (Ka_zz_w1+Ka_zz_w2)];

Ka_disc=double(Ka);
Ks_disc=matlabFunction(Ks);


%% Divergence pressure and velocity

% equivalent straight wing divergence speed
fileID = fopen('straight_div_speed_ms.txt','r');
if fileID == -1
    error('No input straight wing divergence speed!');    % run DivergenceStraight.m before
else
    Ud_straight_ms = fscanf(fileID,'%f');                  % read dynamic divergence pressure in output file
    fclose(fileID);
end

% Initial guess
Cw = 0;
Ud_search_ms = 0;
Ud_tailoring_vect_ms = [];
Cw_increment = 100;

while Ud_search_ms < Ud_straight_ms
    [V, dynPd_Pa] = eig(Ks_disc(Cw), Ka_disc);
    Ud_Tail_ms = sqrt( (diag(dynPd_Pa)*2) / (rho_kgm3) );

    for i=1:(dimN_theta+dimN_z)             % elimination of complex eigenvalues
        if imag(Ud_Tail_ms(i))~=0
            Ud_Tail_ms(i)=inf;
        end
    end

    if min(Ud_Tail_ms) > Ud_straight_ms
        Ud_search_ms = min(Ud_Tail_ms);
    end
    
    Ud_tailoring_vect_ms = [Ud_tailoring_vect_ms min(Ud_Tail_ms)];
    Cw = Cw + Cw_increment;       % iteration
end

Cw_search = Cw-Cw_increment;      % effective Cw value
fprintf('Aeroelastic tailoring positive coupling term: %f Nm^2\n\n',Cw_search);

%% Plot tailoring speed
Cw_vect = (0:Cw_increment:Cw_search);

figure (3)
plot(Cw_vect, Ud_tailoring_vect_ms,LineWidth=1)
hold on
grid on
plot(Cw_vect, Ud_straight_ms*ones(length(Cw_vect),1),LineWidth=1)
hold off
title('Aeroelastic Tailoring')
xlabel('Cw [Nm^2]')
ylabel('U [m/s]')
