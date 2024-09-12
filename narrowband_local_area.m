addpath functions
%%
ch.Ptx = 0.1; % transmit power
ch.wt = 8; % top wall
ch.wb = -8; % bottom wall
ch.rho = 5; % max number of reflections
ch.fc = 3.6e9; % carrier frequency
ch.B = 100e6; % system bandwidth
ch.Q = 1024; % number of frequency bins/OFDM subcarriers
lambda_fc = 3e8/ch.fc; % carrier wavelength
dr = lambda_fc/4; % antenna spacing for MIMO-only (you can change the factor 4)
ch.beta = 2*pi/lambda_fc; % carrier wavenumber
ch.A = ch.Ptx*lambda_fc^2/(4*pi)^2; % pathloss constant
ch.path_params.Npaths = 8; % number of cluster paths per main path
ch.path_params.max_delay = 500e-9; % max path delay, hard cut

%% NarrowBand Local Area
ch.Rx_pos_x = 12:0.01:13; % position of the Client on x
ch.Rx_pos_y = 6:0.01:7; % position of the Client on y

[h_los, h_nlos, theo_K_los, theo_K_nlos, RandTheta] = getNarrowBand(ch); % simulates narrow band mimo channel