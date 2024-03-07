close all
clear

% Eigenmode computation for system:  M*q2 + K*q = 0

run DynamicModel.m

%% Modeshape and Natural_freq computation
[EigVect, EigVal] = eig(Ks_disc, -M_disc);       % eigenmodes and eigenvalues computation
EigVal = diag(EigVal);

EigFreq = (sqrt(EigVal));
EigFreqReal = imag(sqrt(EigVal));               % modes computation

%% Ortogonalization - Gram-schmidt
EigVect_ort = zeros(dimN_theta+dimN_z+dimN_phi);
EigVect_ort(:,1) = EigVect(:,1) / sqrt(M_disc(1,1));            % u1=v1 & normalization on modal mass
j=1;

% Ortogonalization wrt M
for i=(1:(dimN_theta+dimN_z+dimN_phi-1))
    EigVect_ort(:,i+1) = EigVect(:,i+1);    %u2=v2
    for j=(1:i)
        alpha = (EigVect(:,i+1)' * M_disc * EigVect_ort(:,j)) / (EigVect_ort(:,j)' * M_disc * EigVect_ort(:,j));    %proj u1
        EigVect_ort(:,i+1) = EigVect_ort(:,i+1) - ( alpha .* EigVect_ort(:,j) ) ;
    end
end

% diagonalization check
M_modal_mtrx    = EigVect_ort' * M_disc * EigVect_ort;
M_modal         = diag(diag(M_modal_mtrx));        % Modal mass matrix
Ks_modal_mtrx   = EigVect_ort' * Ks_disc * EigVect_ort;
Ks_modal        = diag(diag(Ks_modal_mtrx));       % Modal structural matrix

%% Text output
[EigFreq_sort_sorted, EigFreq_sorted_idx] = sort(EigFreqReal,'ascend');     % eigenfrequencies sorting

% result print in command window
fprintf('First 10 modes of the system are:  ');
fprintf('%f  ',EigFreq_sort_sorted(1:10,1));
fprintf('\n\n');


