close all
clear
clc

%% 1) Modes computation
fprintf('MODES COMPUTATION\n');
run Eigenmode.m                     % swept wing divergence computation script

%% 2) Flutter computation
% WARNING: To run this section must be run part(1) before
fprintf('FLUTTER COMPUTATION\n');
run Flutter.m                       % flutter speed computation
run FlutterPlot.m                   % flutter condition plots


%% 3) Response computation
% WARNING: To run this section must be run part(2) before
fprintf('FREQUENCY RESPONSE COMPUTATION\n');
run FreqResponse.m                  % frequency response computation
