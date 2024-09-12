addpath functions
clear; close all; clc;
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
x_res = 0.01;
y_res = 0.01;
ch.Rx_pos_x = 12:x_res:13; % position of the Client on x
ch.Rx_pos_y = 6:y_res:7; % position of the Client on y
dimx = length(ch.Rx_pos_x);
dimy = length(ch.Rx_pos_y);
local_area_x = mean(ch.Rx_pos_x);
local_area_y = mean(ch.Rx_pos_y);
los = 1;

%% Plot the mean |h| for different subcarriers
band_gap = ch.B/(ch.Q-1);
plot_band_gap_ratio = 100;
subcarriers_to_plot = ch.fc-ch.B/2:plot_band_gap_ratio*band_gap:ch.fc+ch.B/2;
am_subcarriers = length(subcarriers_to_plot);
am_realizations = 50;
power_matrix_los = zeros(am_subcarriers, am_realizations);
power_matrix_nlos = zeros(am_subcarriers, am_realizations);
K_matrix_los = zeros(am_subcarriers, am_realizations);
K_matrix_nlos = zeros(am_subcarriers, am_realizations);
for subcarrier_index = 1:am_subcarriers
    ch.fc = subcarriers_to_plot(subcarrier_index);
    for realization_index = 1:am_realizations
        [h_los, h_nlos, K_los, K_nlos,~] = getNarrowBand(ch); % simulates narrow band mimo channel
        h_los_abs_mean = mean(abs(h_los),"all");
        h_nlos_abs_mean = mean(abs(h_nlos),"all");
        K_los_abs_mean = mean(abs(K_los),"all");
        K_nlos_abs_mean = mean(abs(K_nlos),"all");
        power_matrix_los(subcarrier_index, realization_index) = h_los_abs_mean;
        power_matrix_nlos(subcarrier_index, realization_index) = h_nlos_abs_mean;
        K_matrix_los(subcarrier_index, realization_index) = K_los_abs_mean;
        K_matrix_nlos(subcarrier_index, realization_index) = K_nlos_abs_mean;
    end
end
power_array_los = mean(power_matrix_los, 2);
power_array_nlos = mean(power_matrix_nlos,2);
K_array_los = mean(K_matrix_los,2);
K_array_nlos = mean(K_matrix_nlos,2);

%finally, some hard values
power_mean_los = mean(power_array_los);
power_mean_nlos = mean(power_array_nlos);
K_mean_los = mean(K_array_los);
K_mean_nlos = mean(K_array_nlos);

save("narrowband_2_1_partone.mat");

f = figure;
title("average |h| at different subcarriers", "\Deltaf = " + num2str((subcarriers_to_plot(2)-subcarriers_to_plot(1))/1e6,3) + " MHz = " + num2str(plot_band_gap_ratio) + "* band gap")
xlabel("f [Hz]")
ylabel('|h|')
hold on
plot(subcarriers_to_plot, power_array_los, '-o', 'DisplayName', 'los');
plot(subcarriers_to_plot, power_array_nlos, '-o', 'DisplayName', 'nlos');
legend('show');
exportgraphics(f,'Report/Figures/2_1/CTF_fc.png', 'Resolution', 1000);

f_2 = figure;
title("average theoretical Rice factor K at different subcarriers", "\Deltaf = " + num2str((subcarriers_to_plot(2)-subcarriers_to_plot(1))/1e6,3) + " MHz = " + num2str(plot_band_gap_ratio) + "* band gap")
xlabel("f [Hz]")
ylabel('K')
hold on
plot(subcarriers_to_plot, K_array_los, '-o', 'DisplayName', 'los');
plot(subcarriers_to_plot, K_array_nlos, '-o', 'DisplayName', 'nlos');
legend('show');
exportgraphics(f_2,'Report/Figures/2_1/K_fc.png', 'Resolution', 1000);

%% obtain for rest of 2.1
ch.fc = 3.6e9; % restore
[h_los, h_nlos, theo_K_los, theo_K_nlos, ~] = getNarrowBand(ch); % simulates narrow band mimo channel
%h : dimy*dimx
%----> x
%|
%|
%|
%v
% y

%RandTheta: 101x101x11x2: pos_y*pos_x*?*?

K_theo_mean_los = mean(mean(theo_K_los));
K_theo_mean_nlos = mean(mean(theo_K_nlos));

%% obtain PDF and CDF

h_los_reshaped = reshape(h_los,dimx*dimy,1); %histogram functions take vectors
h_nlos_reshaped = reshape(h_nlos,dimx*dimy,1); %histogram functions take vectors

figure
t_1 = tiledlayout(2,1);
title(t_1, "histograms of |h| over local area", "around (x,y) = (" + num2str(local_area_x) + "," + num2str(local_area_y) + ")")
nexttile(t_1)
histfit(abs(h_los_reshaped),ch.Q,'rician');
title("histogram of |h_{los}| w fitted rician pdf")
xlabel("|h|")
ylabel("#occurrences in local area")
nexttile(t_1)
histfit(abs(h_nlos_reshaped),ch.Q,'rayleigh');
title("histogram of |h_{nlos}| w fitted rayleigh pdf")
xlabel("|h|")
ylabel("#occurrences in local area")
hold off
exportgraphics(t_1,'Report/Figures/2_1/histograms.png', 'Resolution', 1000);


data_los = fitdist(abs(h_los_reshaped),'rician'); %look into for first K-value
data_nlos = fitdist(abs(h_nlos_reshaped),'rayleigh');
K_1_los = data_los.s/data_los.sigma;
figure
bin_resolution = 10000;
hist_los = histogram(abs(h_los_reshaped), bin_resolution); %we need this struct for the x-axis (code: https://nl.mathworks.com/matlabcentral/answers/453732-changing-histogram-to-pdf)
bin_centers_los = hist_los.BinEdges(1:end-1) + (hist_los.BinWidth/2); %code: https://nl.mathworks.com/matlabcentral/answers/453732-changing-histogram-to-pdf
hist_nlos = histogram(abs(h_nlos_reshaped), bin_resolution); %we need this struct for the x-axis (code: https://nl.mathworks.com/matlabcentral/answers/453732-changing-histogram-to-pdf)
bin_centers_nlos = hist_nlos.BinEdges(1:end-1) + (hist_nlos.BinWidth/2); %code: https://nl.mathworks.com/matlabcentral/answers/453732-changing-histogram-to-pdf
p_los = pdf(data_los,bin_centers_los);
p_norm_los = p_los/sum(p_los);
p_nlos = pdf(data_nlos,bin_centers_los);
p_norm_nlos = p_nlos/sum(p_nlos);
%idea: cdf cumulative result of pdf
c_los = cumsum(p_norm_los);
c_nlos = cumsum(p_norm_nlos);
figure
t_2 = tiledlayout(2,1);
title(t_2, "propability and cumulative density functions of |h| over local area", "around (x,y) = (" + num2str(local_area_x) + "," + num2str(local_area_y) + ")")
nexttile(t_2)
hold on
plot(bin_centers_los,1000*p_norm_los, 'DisplayName', 'pdf*1000' )
plot(bin_centers_los,c_los, 'DisplayName', 'cdf')
title("Empirical result of PDF and CDF in Local area, LOS")
xlabel("|h|")
legend('show')
nexttile(t_2)
hold on
plot(bin_centers_nlos,1000*p_norm_nlos, 'DisplayName', 'pdf*1000' )
plot(bin_centers_nlos,c_nlos, 'DisplayName', 'cdf')
title("Empirical result of PDF and CDF in Local area, NLOS")
xlabel("|h|")
legend('show')
exportgraphics(t_2,'Report/Figures/2_1/pdf_cdf.png', 'Resolution', 1000);

%% spatial correlation along x-axis (results in R-array ifo x-position)
%----> x
%|
%|
%|
%v
% y

n = size(h_los,2); %amount of h_values along x-axis
R_x_los = zeros(n,2*dimy-1);
R_x_nlos = zeros(n,2*dimy-1);

for y_index = 1:dimy
    R_x_los(y_index,:) = xcorr(h_los(y_index,:));
    R_x_nlos(y_index,:) = xcorr(h_nlos(y_index,:));
end
R_x_los = R_x_los(:,dimy:end);
%R_x_los = R_x_los(1,:);
R_x_los = mean(R_x_los,1);
R_x_los= R_x_los./ max(abs(R_x_los)); %normalisation

R_x_nlos = R_x_nlos(:,dimy:end);
%R_x_nlos = R_x_nlos(1,:);
R_x_nlos = mean(R_x_nlos,1);
R_x_nlos= R_x_nlos./ max(abs(R_x_nlos)); %normalisation

threshold = 0.7;
horizontal_curve = threshold * ones(size(R_x_los));
x_axis = ch.Rx_pos_x - ch.Rx_pos_x(1);
f_Rx = figure;
hold on
title("Narrowband Channel spatial correlation, x-direction", "around ("+num2str(local_area_x)+","+num2str(local_area_y)+")")

xlabel("\Deltax [m]")
ylabel("R(\Deltax) (mean along y-axis)")
plot(x_axis,abs(R_x_los),'DisplayName',strcat('R(\Deltax)_{LOS}'), 'Color', 'r')
plot(x_axis,abs(R_x_nlos),'DisplayName',strcat('R(\Deltax)_{NLOS}'), 'Color', 'g')
plot(x_axis,horizontal_curve,'DisplayName',strcat('threshold'))
point_los = get_threshold_value_index(R_x_los, threshold);
point_nlos = get_threshold_value_index(R_x_nlos, threshold);
plot((point_los-1)*(x_axis(2)-x_axis(1)), abs(R_x_los(point_los)), 'xr', 'DisplayName',strcat('\Deltax =', num2str(x_axis(point_los)),"m"), 'Color','r')
plot((point_nlos-1)*(x_axis(2)-x_axis(1)), abs(R_x_nlos(point_nlos)), 'xr', 'DisplayName',strcat('\Deltax =', num2str(x_axis(point_nlos)),"m"), 'Color', 'g')
xlim([x_axis(1) x_axis(end)]);
legend('show')
hold off
exportgraphics(f_Rx,'Report/Figures/2_1/Rx.png', 'Resolution', 1000);

%% spatial correlation along y-axis

%----> x
%|
%|
%|
%v
% y

n = size(h_los,1); %amount of h_values along the y-axis
R_y_los = zeros(n,2*dimx-1);
R_y_nlos = zeros(n,2*dimx-1);

for x_index = 1:dimx
    R_y_los(x_index,:) = xcorr(h_los(:,x_index));
    R_y_nlos(x_index,:) = xcorr(h_nlos(:,x_index));
end
R_y_los = R_y_los(:,dimy:end);
%R_y = R_y(1,:);
R_y_los = mean(R_y_los,1);
R_y_los= R_y_los./ max(abs(R_y_los)); %normalisation

R_y_nlos = R_y_nlos(:,dimy:end);
%R_y = R_y(1,:);
R_y_nlos = mean(R_y_nlos,1);
R_y_nlos= R_y_nlos./ max(abs(R_y_nlos)); %normalisation

threshold = 0.7;
horizontal_curve = threshold * ones(size(R_y_los));
y_axis = ch.Rx_pos_y - ch.Rx_pos_y(1);
f_Ry = figure;
hold on
title("Narrowband Channel spatial correlation, y-direction", "around ("+num2str(local_area_x)+","+num2str(local_area_y)+")")

xlabel("\Deltay [m]")
ylabel("R(\Deltay) (mean along x-axis)")
plot(y_axis,abs(R_y_los),'DisplayName',strcat('R(\Deltay)_{LOS}'), 'Color', 'r')
plot(y_axis,abs(R_y_nlos),'DisplayName',strcat('R(\Deltay)_{NLOS}'), 'Color', 'g')
plot(y_axis,horizontal_curve,'DisplayName',strcat('threshold'))
point_los = get_threshold_value_index(R_y_los, threshold);
point_nlos = get_threshold_value_index(R_y_nlos, threshold);
plot((point_los-1)*(y_axis(2)-y_axis(1)), abs(R_y_los(point_los)), 'xr', 'DisplayName',strcat('\Deltay =', num2str(y_axis(point_los)),"m"), 'Color','r')
plot((point_nlos-1)*(y_axis(2)-y_axis(1)), abs(R_y_nlos(point_nlos)), 'xr', 'DisplayName',strcat('\Deltay =', num2str(y_axis(point_nlos)),"m"), 'Color', 'g')
xlim([y_axis(1) y_axis(end)]);
legend('show')
hold off
exportgraphics(f_Ry,'Report/Figures/2_1/Ry.png', 'Resolution', 1000);

%% PAS
theta_axis = linspace(0,pi,1000);
u_axis = ch.beta*cos(theta_axis);
S_x_los = zeros(size(u_axis));
S_x_nlos = zeros(size(u_axis));
for theta_index = 1:length(theta_axis)
    u_curr = u_axis(theta_index);
    S_curr = 0;
    for R_index = 1:length(R_x_los)
        R_curr = R_x_los(R_index);
        S_curr = S_curr + R_curr*exp(-1j*u_curr*x_res*(R_index-1));
    end
    S_curr = abs(S_curr);
    S_x_los(theta_index) = S_curr;

    S_curr = 0;
    for R_index = 1:length(R_x_nlos)
        R_curr = R_x_nlos(R_index);
        S_curr = S_curr + R_curr*exp(-1j*u_curr*x_res*(R_index-1));
    end
    S_curr = abs(S_curr);
    S_x_nlos(theta_index) = S_curr;
end

 % Total incident power
P_x_los = trapz(u_axis, abs(S_x_los));
P_x_nlos = trapz(u_axis, abs(S_x_nlos));
% Mean AoA
u_x_m_los = trapz(u_axis, u_axis.*abs(S_x_los))/P_x_los;
u_x_m_nlos = trapz(u_axis, u_axis.*abs(S_x_nlos))/P_x_nlos;
% Angular spread
sigma_x_u_los = sqrt(trapz(u_axis, u_axis.^2 .* abs(S_x_los))/P_x_los - u_x_m_los^2);
sigma_x_u_nlos = sqrt(trapz(u_axis, u_axis.^2 .* abs(S_x_nlos))/P_x_nlos - u_x_m_nlos^2);

S_y_los = zeros(size(u_axis));
S_y_nlos = zeros(size(u_axis));
for theta_index = 1:length(theta_axis)
    u_curr = u_axis(theta_index);
    S_curr = 0;
    for R_index = 1:length(R_y_los)
        R_curr = R_y_los(R_index);
        S_curr = S_curr + R_curr*exp(-1j*u_curr*y_res*(R_index-1));
    end
    S_curr = abs(S_curr);
    S_y_los(theta_index) = S_curr;

    S_curr = 0;
    for R_index = 1:length(R_y_nlos)
        R_curr = R_y_nlos(R_index);
        S_curr = S_curr + R_curr*exp(-1j*u_curr*y_res*(R_index-1));
    end
    S_curr = abs(S_curr);
    S_y_nlos(theta_index) = S_curr;
end

 % Total incident power
P_y_los = trapz(u_axis, abs(S_y_los));
% Mean AoA
u_y_m_los = trapz(u_axis, u_axis.*abs(S_y_los))/P_y_los;
% Angular spread
sigma_y_u_los = sqrt(trapz(u_axis, u_axis.^2 .* abs(S_y_los))/P_y_los - u_y_m_los^2);

 % Total incident power
P_y_nlos = trapz(u_axis, abs(S_y_nlos));
% Mean AoA
u_y_m_nlos = trapz(u_axis, u_axis.*abs(S_y_nlos))/P_y_nlos;
% Angular spread
sigma_y_u_nlos = sqrt(trapz(u_axis, u_axis.^2 .* abs(S_y_nlos))/P_y_nlos - u_y_m_nlos^2);

f_PAS = figure;
hold on
plot(theta_axis, 10*log10(abs(S_x_los)), 'DisplayName',"S_{x,LOS}")
plot(theta_axis, 10*log10(abs(S_y_los)), 'DisplayName', "S_{y,LOS}")
plot(theta_axis, 10*log10(abs(S_x_nlos)), 'DisplayName',"S_{x,NLOS}")
plot(theta_axis, 10*log10(abs(S_y_nlos)), 'DisplayName', "S_{y,NLOS}")

title("PAS along x- and y-axis", "around ("+num2str(local_area_x)+","+num2str(local_area_y)+")")
xlabel("\theta [radians]")
ylabel("10log_{10}(PAS) [dB]")
legend("show")
exportgraphics(f_PAS,'Report/Figures/2_1/PAS.png', 'Resolution', 1000)
%% get to MIMOfunction
ch.Rx_pos_x = 12; % position of the Client on x
ch.Rx_pos_y = 6; % position of the Client on y
ch.Ntx = 4; % number of antennas at the Base Station
ch.Nrx = 4; % number of antennas at the Client
ch.RxArrayOrient = 'x'; % 'x' or 'y' % orientation of the linear array
ch.AntSpacing = dr;

am_antennae = 1:4;
eigenvalues_los_store = cell(1,am_antennae(end)); % amount of eigenvalues will vary!
eigenvalues_nlos_store = cell(1,am_antennae(end)); % amount of eigenvalues will vary!

for am_antenna = am_antennae
    ch.Ntx = am_antenna; % number of antennas at the Base Station
    ch.Nrx = am_antenna; % number of antennas at the Client
    [h_los, h_nlos, ~] = getNarrowBandMIMO(ch); % simulates narrow band mimo channel
    [U_los,S_los,V_los] = svd(h_los);
    [U_nlos, S_nlos, V_nlos] = svd(h_nlos);
    S_los_diag = zeros(1,min(size(S_los)));
    for i = 1:min(size(S_los))
        S_los_diag(i) = S_los(i,i);
    end
    eigenvalues_los_store{am_antenna} = S_los_diag;

    S_nlos_diag = zeros(1,min(size(S_nlos)));
    for i = 1:min(size(S_nlos))
        S_nlos_diag(i) = S_nlos(i,i);
    end
    eigenvalues_nlos_store{am_antenna} = S_nlos_diag;
end

%% Plot eigenvalues
colours = ['r','b','g','k'];

figure
t_2 = tiledlayout(1,2);
title(t_2, "CTF eigenvalues i.f.o #antennae", "around coord (" + num2str(ch.Rx_pos_x) + "," + num2str(ch.Rx_pos_y) + ")")
nexttile
hold on
for am_antenna = am_antennae
    current_eigenvalues_los = eigenvalues_los_store{am_antenna};
    current_colour = colours(am_antenna);
    for current_eigenvalue_los = current_eigenvalues_los
        plot(am_antenna, current_eigenvalue_los, 'o', 'Color', current_colour)
    end
end
title("LOS")
xlabel("n° antennae Tx and Rx")
ylabel("CTF eigenvalues")

nexttile
hold on
for am_antenna = am_antennae
    current_eigenvalues_nlos = eigenvalues_nlos_store{am_antenna};
    current_colour = colours(am_antenna);
    for current_eigenvalue_nlos = current_eigenvalues_nlos
        plot(am_antenna, current_eigenvalue_nlos, 'o', 'Color', current_colour)
    end
end
title("NLOS")
xlabel("n° antennae Tx and Rx")
ylabel("CTF eigenvalues")

%% save for later
save("2.1_workbench.mat")

%% funcs
function point = get_threshold_value_index(array, threshold)
    for i = 1:length(array)
        if (abs(array(i)) < threshold)
            point = i-1;
            return
        end
    end
    point = length(array);
    return
end

function P = calculate_incident_power(u, S)
    P = trapz(u,S);
end

function u_m = calculate_mean_AoA(u, S, P)
    Y = u.*S;
    u_m = trapz(u, Y) / P;
end

function sigma_u = calculate_angular_spread(u,S,P,u_m)
    Y = (u.^2).*S;
    sigma_u = abs(trapz(u, Y))/P;
    sigma_u = sigma_u - u_m^2;
    sigma_u = sqrt(sigma_u);
end

function coordinates = calculate_coordinates_around_center(cent_x,cent_y,am_x,am_y,dist_x,dist_y)
    coordinates = zeros(am_x*am_y, 2);
    if mod(am_x,2) == 0 %even amount
        x_range = cent_x-((am_x)/2-1)*dist_x:dist_x:cent_x+(am_x/2)*dist_x;
    else
        x_range = cent_x-((am_x-1)/2)*dist_x:dist_x:cent_x+((am_x-1)/2)*dist_x;
    end

    if mod(am_y,2) == 0 %even amount
        y_range = cent_y-(am_y/2-1)*dist_y:dist_y:cent_y+(am_y/2)*dist_y;
    else
        y_range = cent_y-((am_y-1)/2)*dist_y:dist_y:cent_y+((am_y-1)/2)*dist_y;
    end
    index = 1;
    for x = x_range
        for y = y_range
            coordinates(index, :) = [x, y];
            index = index + 1;
        end
    end
end