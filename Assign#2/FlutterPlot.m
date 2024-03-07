% Flutter plots for system:  M*q'' + C*q' + K*q = 0
% WARNING: To run this script must be previous computed 'Flutter.m'

%% Root locus

figure(1)
plot(real(Eigenvalues_vect), imag(Eigenvalues_vect), LineStyle="none",Marker="o")
grid on
title ('Locus root')
ylabel ('immaginary')
xlabel ('real')

%% G-damping

%Plot for 11 Modeshapes
figure(2)
plot(U_inf_ms, g_damp(1,:),Marker="x",SeriesIndex=1)
hold on
plot(U_inf_ms, g_damp(3,:),Marker="x",SeriesIndex=2)
plot(U_inf_ms, g_damp(7,:),Marker="x",SeriesIndex=3)
plot(U_inf_ms, g_damp(9,:),Marker="x",SeriesIndex=4)
plot(U_inf_ms, g_damp(11,:),Marker="x",SeriesIndex=5)
plot(U_inf_ms, g_damp(13,:),Marker="x",SeriesIndex=6)
plot(U_inf_ms, g_damp(15,:),Marker="x",SeriesIndex=7)
plot(U_inf_ms, g_damp(17,:),Marker="x",SeriesIndex=8)
plot(U_inf_ms, g_damp(19,:),Marker="x",SeriesIndex=9)
plot(U_inf_ms, g_damp(21,:),Marker="x",SeriesIndex=10)
plot(U_inf_ms, g_damp(23,:),Marker="x",SeriesIndex=11)
grid on
axis on
title ('g-damping plot')
ylabel ('g-damp')
xlabel ('U m/s')
hold off