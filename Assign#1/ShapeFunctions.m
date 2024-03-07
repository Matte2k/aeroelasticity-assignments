% WARNING - must be previously defined the symbolic variable 'y_mark_m' and 'Data.m'

%% MultipleSHAPE function definition
     
% Theta shapes vector definition
num_Ntheta = 2;
N_theta =   [ cos([1 : num_Ntheta].*y_mark_m.*(pi/bw_mark_m))-1, ...
                sin( (pi/(2*bw_mark_m) + [0 : num_Ntheta-1].*(pi/bw_mark_m)).*y_mark_m ) ];

% Z shapes vector definition
num_Nz = 2;
N_z =       [ cos( (pi/(2*bw_mark_m) + [0 : num_Nz - 1].*(pi/bw_mark_m))*y_mark_m )-1, ...
                sin( [1 : num_Nz] .* y_mark_m *(pi/bw_mark_m) ) - ([1 : num_Nz].*y_mark_m*(pi/bw_mark_m)) ];


dimN_theta = length(N_theta);       % Theta shapes vector dimension
dimN_z = length(N_z);               % Z shapes vector dimension

% derivative of shapes funciton vectors
dN_theta = diff(N_theta);
dN_z = diff(N_z);
d2N_z = diff(N_z,2);

