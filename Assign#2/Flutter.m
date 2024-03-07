% Flutter computation for system:  M*q'' + C*q' + K*q = 0

%% Speed vector definition

U_inf_ms = 0:1:150;                         % speed vector [ms]
eps_flutter = 1e-6;                         % limit over flutter


%% Modal semplification of the system

Minv_modal = inv(M_modal);      % inverse of modal mass
I_modal = eye(Sys_dim);         % identity matrix for state space system

Cs_modal = diag(2 .* diag(M_modal) .* damp_fraction .* EigFreqReal);    % modal damping matrix  

Eigenvalues_vect = zeros(2*Sys_dim,(length(U_inf_ms)));                   % eigenvalues vector initialization

for i = 1:(length(U_inf_ms))
    Ca_modal = EigVect_ort' * Ca_disc * EigVect_ort;                    % modal aerodynamic damping matrix
    C_modal = (Cs_modal) - ( 1/2 * rho_kgm3 * U_inf_ms(i)) * (Ca_modal)  ;       % modal total damping
    K_modal = EigVect_ort' * (Ks_disc - (1/2 * rho_kgm3 * U_inf_ms(i)^2) * (Ka_disc)) * EigVect_ort;    % modal total stiffness

    % modal state space system definition
    StateSpaceSys = [   -Minv_modal*C_modal,    -Minv_modal*K_modal; ...
                        I_modal,                zeros(Sys_dim)  ];

    Sys_Eig = eig(StateSpaceSys);                   % state space eigenvalues
    Eigenvalues_vect(:,i) = Sys_Eig;                  % eigenvalues vector for each speed
end

g_damp = real(Eigenvalues_vect)./imag(Eigenvalues_vect);    % g-damping computation

% Im==0 modes elimination
for j = 1:(height(g_damp))
    if norm(real(Eigenvalues_vect(j,1))) < (eps_flutter)
        g_damp(j,:) = NaN;
    end
end

%% Flutter speed computation 

U_flutter_ms = zeros(Sys_dim*2,1);
U_flutter_ms(2:2:22) = inf;

for i=1:2:(2*Sys_dim)
    for j=1:length(U_inf_ms)
        if imag(Eigenvalues_vect(i,j))~=0 && g_damp(i,j)>eps_flutter && U_flutter_ms(i)==0
            U_flutter_ms(i) = U_inf_ms(j);
        end
    end
end


%% Text output
[U_first_flutter_ms, flutter_mode_idx] = min(U_flutter_ms(U_flutter_ms>0));     % minimum flutter speed>0 selection 

% result print in command window
fprintf('First flutter speed is: %f m/s \n\n', U_first_flutter_ms);
