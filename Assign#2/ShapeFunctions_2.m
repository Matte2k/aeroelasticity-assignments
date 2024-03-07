% WARNING - must be previously defined the symbolic variable 'y_mark_m', 'x_m' and 'Data2.m'

%% MultipleSHAPE function definition
num_Ntheta = 2;
num_Nz = 2;
num_Nphi = 1;

%% rigid in
% Wing-bend shape function (with rigid)
N_z =       [1, cos( (pi/(2*bw_mark_m) + [0 : num_Nz - 1].*(pi/bw_mark_m))*y_mark_m )-1,    sin( [1 : num_Nz] .* y_mark_m *(pi/bw_mark_m) ) - ([1 : num_Nz].*y_mark_m*(pi/bw_mark_m)) ];     %INTERO AEREO condizioni BC derivata seconda

% Wing-twist shape function (with rigid)
N_theta =   [1, cos([1 : num_Ntheta].*y_mark_m.*(pi/bw_mark_m))-1,        sin( (pi/(2*bw_mark_m) + [0 : num_Ntheta-1].*(pi/bw_mark_m)).*y_mark_m ) ];                      %INTERO AEREO condizioni

% Fuselage-twist shape function (with rigid)
N_phi =   [1, 1-x_m/norm(PTcstrC_m)];

dimN_theta = length(N_theta);
dimN_z = length(N_z);
dimN_phi = length(N_phi);
Sys_dim = dimN_z + dimN_theta + dimN_phi;

% Derivative
dN_theta = diff(N_theta);
dN_z = diff(N_z);
d2N_z = diff(N_z,2);
dN_phi = diff(N_phi);

