% Root moment computed with symbolic computation
% WARNING: To run this script must be previous computed 'Rollrate.m'

%% ALEIRON moment at the root 

AleironMoment_w1_Nm =   + int( q_roll_Pa .* cw_m .* CLbetaA .* Beta_rad .* cos(SweepAng_rad_w1).^2 .* y_mark_m, y_mark_m, ba_mark_m, bw_mark_m);
AleironMoment_w2_Nm =   + int( q_roll_Pa .* cw_m .* CLbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2).^2 .* y_mark_m, y_mark_m, -bw_mark_m, -ba_mark_m);
AleironMoment_Nm = double(AleironMoment_w1_Nm+AleironMoment_w2_Nm);
fprintf('Aleiron moment at wingroot: %f Nm\n',AleironMoment_Nm);


%% WING moment at the root

WingMoment_w1_Nm =    + int ( theta_mark .* q_roll_Pa .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w1).^3, y_mark_m, 0, bw_mark_m) ...
                      - int ( dz .* q_roll_Pa .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w1).^2 .* sin(SweepAng_rad_w1) , y_mark_m, 0, bw_mark_m) ...                                                            
                      - int ( p_rads .* q_roll_Pa .* cw_m .* CLalfaW .* (y_mark_m).^2 .* (1 / Uinf_ms) .* cos(SweepAng_rad_w1).^3 , y_mark_m, 0, bw_mark_m) ...                                        
                      + int ( q_roll_Pa .* cw_m .* CLbetaA .* Beta_rad .* cos(SweepAng_rad_w1).^2 .* y_mark_m, y_mark_m, ba_mark_m, bw_mark_m);
WingMoment_w2_Nm =    + int ( theta_mark .* q_roll_Pa .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w2).^3, y_mark_m, -bw_mark_m, 0) ...
                      - int ( dz .* q_roll_Pa .* cw_m .* CLalfaW .* y_mark_m .* cos(SweepAng_rad_w2).^2 .* sin(SweepAng_rad_w2) , y_mark_m, -bw_mark_m, 0) ...                                                            
                      - int ( p_rads .* q_roll_Pa .* cw_m .* CLalfaW .* (y_mark_m).^2 .* (1 / Uinf_ms) .* cos(SweepAng_rad_w2).^3 , y_mark_m, -bw_mark_m, 0) ...                                        
                      + int ( q_roll_Pa .* cw_m .* CLbetaA .* (-Beta_rad) .* cos(SweepAng_rad_w2).^2 .* y_mark_m, y_mark_m, -bw_mark_m, -ba_mark_m);
WingMoment_Nm = double(WingMoment_w1_Nm + WingMoment_w2_Nm);
fprintf('Aerodynamic moment at wingroot: %f Nm\n',WingMoment_Nm);


%% CANARD moment at the root

CanardMoment_w1_Nm = - int ( p_rads .* q_roll_Pa .* cc_m .* CLalfaC .* (y_mark_m).^2 .* (1 / Uinf_ms), y_mark_m, 0, bc_m);
CanardMoment_w2_Nm = - int ( p_rads .* q_roll_Pa .* cc_m .* CLalfaC .* (y_mark_m).^2 .* (1 / Uinf_ms), y_mark_m, -bc_m, 0);
CanardMoment_Nm = double(CanardMoment_w1_Nm+CanardMoment_w2_Nm);
fprintf('Aerodynamic moment at canardroot: %f Nm\n\n',CanardMoment_Nm);


%% FUSELAGE moment plot

x_plot1 = [-lf_m : 0.01 : PTcstrC_m(1,1)];
TorsionMoment_vect1 = zeros(length(x_plot1),1);

x_plot2 = [PTcstrC_m(1,1) : 0.01 : 0];
TorsionMoment_vect2 = CanardMoment_Nm .* ones(length(x_plot2),1);

figure(2)
hold on
plot([x_plot1 x_plot2],zeros(1,length([x_plot1 x_plot2])),'k',LineWidth=3)      % fuselage plot
plot(x_plot1, TorsionMoment_vect1,'r',LineWidth=1)
plot([1 1]*x_plot1(length(x_plot1)),[0,TorsionMoment_vect2(1)],'r',LineWidth=1)
plot([0 0],[0,TorsionMoment_vect2(1)],'r',LineWidth=1)
plot(x_plot2, TorsionMoment_vect2,'r',LineWidth=1)
xlim([-10 1])
ylim([-600 100])
title('Torsional moment plot')
xlabel('x')
ylabel('M [Nm]')
legend('Fuselage', 'Torsional moment',Location='southwest')
grid on
hold off

