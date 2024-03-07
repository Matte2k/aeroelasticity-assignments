% Displacements plot in elastic roll condition

% Wing twist displacement function
theta_mark = N_theta * res([1:dimN_theta],1);           % theta symbolic function definition
theta_mark_fun = matlabFunction(theta_mark);            % theta Matlab function conversion

% Wing twist displacement function
z = N_z * res([(dimN_theta+1):(dimN_z+dimN_theta)],1);  % z symbolic function definition
dz= diff(z,1);                                          % z derivative use later in 'RootMoment.m'
z_fun = matlabFunction(z);                              % z Matlab function conversion

y_plot = [-bw_mark_m:0.01:bw_mark_m];      % y axis definition
theta_scale = 1;                           % scale factor to visualize better theta displacement

figure (1)
plot(y_plot,(theta_scale.*theta_mark_fun(y_plot)),'LineWidth',1)    % theta function plot
grid on
hold on
plot(y_plot,z_fun(y_plot),'LineWidth',1)                            % z function plot
legend('theta', 'z')
title ('Wing Displacements in roll condition')
drawnow