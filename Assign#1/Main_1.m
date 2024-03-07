close all
clear
clc

%% 1) Divergence speed computation
fprintf('DIVERGENCE COMPUTATION\n');
run Divergence.m                    % swept wing divergence computation script

%% 2) Elastic and Rigid Rollrate
fprintf('ROLLRATE COMPUTATION\n');
run Rollrate_Rigid.m                % rigid roll computation script
run Rollrate.m                      % elastic roll computation script
run Roll_Plot_displacements.m       % displacement in roll plot script

%% 3-4) Torsional moment diagram
% WARNING: To run this section must be run part(2) before
fprintf('MOMENT at WINGROOT COMPUTATION\n');
run RootMoment.m                    % moments computation and plot script

%% 5) Aeroelastic tailoring
fprintf('AEROELASTIC TAILORING COMPUTATION\n');
run DivergenceStraight.m            % straight wing divergence computation script
run AeroelasticTailoring.m          % aeroelastic tailorin script