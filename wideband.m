clc; close all; clear all;
addpath functions
%%
ch.Ptx = 0.1;                         % transmit power
ch.wt = 8;                            % top wall
ch.wb = -8;                           % bottom wall
ch.rho = 3;                           % max number of reflections
ch.fc = 3.6e9;                        % carrier frequency
ch.B = 100e6;                         % system bandwidth
ch.Q = 1024;                          % number of frequency bins/OFDM subcarriers
lambda_fc = 3e8/ch.fc;                % carrier wavelength
dr = lambda_fc/4;                     % antenna spacing for MIMO-only (you can change the factor 4)
ch.beta = 2*pi/lambda_fc;             % carrier wavenumber
ch.A = ch.Ptx*lambda_fc^2/(4*pi)^2;   % pathloss constant
ch.path_params.Npaths = 8;            % number of cluster paths per main path (S-V Channel Model)
ch.path_params.max_delay = 500e-9;    % max path delay, hard cut
%% WideBand
ch.Rx_pos_x = 12;
ch.Rx_pos_y = 6;
% subcarrier_spacing= ch.B/ch.Q;       = 97.65625 kHz
% minimum sampling rate required in this case would be 2*97.65625 kHz = 195.3125 kHz.


[H_los, H_nlos, RandTheta] = getWideBand(ch); % simulates wideband siso channel

client_position = 1:5:21;
ctf_am_reflection = [0, 1, 2, 10];
w_freq = 0:ch.B/ch.Q:(ch.B-ch.B/ch.Q); % Frequency

figure

line_colors = {'b', 'r', 'g', 'm'}; % Define line colors

for am_ref_ind = 1:length(ctf_am_reflection)
    all_reff_num = ctf_am_reflection(am_ref_ind);
    ch.rho = all_reff_num;

    subplot(4, 2, am_ref_ind*2 - 1)
    hold on
    % title("CFR for max LOS " + ch.rho + " reflections")
    title("CFR LOS " + ch.rho + " reflection(s)")
    xlabel("Frequency")
    ylabel("Reflections")

    for i = 1:length(client_position)
        distance = client_position(i);
        [~, color_idx] = ismember(distance, [1, 6, 11, 16, 21]);
        
        if color_idx > 0 && color_idx <= length(line_colors)
            line_color = line_colors{color_idx};
        else
            line_color = 'k'; % Use black color for other distances
        end
        
        ch.Rx_pos_x = distance;
        [H_los, ~, ~] = getWideBand(ch);
        plot(w_freq, 10*log10(abs(H_los)), 'DisplayName', strcat('x=', num2str(distance)), 'Color', line_color)
    end
    legend('show')
    hold off

    subplot(4, 2, am_ref_ind*2)
    hold on
    % title("NLOS max reflection " + ch.rho + " reflections")
    title("CFR NLOS " + ch.rho + " reflection(s)")
    xlabel("Frequency")
    ylabel("Reflections")

    for i = 1:length(client_position)
        distance = client_position(i);
        [~, color_idx] = ismember(distance, [1, 6, 11, 16, 21]);
        
        if color_idx > 0 && color_idx <= length(line_colors)
            line_color = line_colors{color_idx};
        else
            line_color = 'k'; % Use black color for other distances
        end
        
        ch.Rx_pos_x = distance;
        [~, H_nlos, ~] = getWideBand(ch);
        plot(w_freq, 10*log10(abs(H_nlos)), 'DisplayName', strcat('x=', num2str(distance)), 'Color', line_color)
    end
    legend('show')
    hold off
end
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ch.Rx_pos_x = 12;                % back to default value
all_reff_num = [0, 5, 10];
n_datapoint = length(H_los);      % amount of datapoints
hard_threshold = 0.7;             % as long as the channel frequency correlation is larger than this value, we deem it relevant
% 9.164634384690092e-08
horizontal_threshold_fixed = hard_threshold * ones(1, floor(n_datapoint/2));
averaging_amount = 150;

% Realization of the CFR of the LOS and NLOS
H_los_buffer = zeros(length(all_reff_num), n_datapoint);
for ctf_am_reflection = 1:length(all_reff_num)
    ch.rho = all_reff_num(ctf_am_reflection);
    
    % averaging H_los, H_nlos because we have a stochastic process
    Total_H_los = 0;
    Total_nH_los = 0;
    for avg_ind = 1:averaging_amount
        [H_los, H_nlos, ~] = getWideBand(ch);
        Total_H_los = Total_H_los + H_los;
        Total_nH_los = Total_nH_los + H_nlos;
    end
    Total_H_los = Total_H_los / averaging_amount;
    Total_nH_los = Total_nH_los / averaging_amount;
    H_los_buffer(ctf_am_reflection, :) = Total_H_los;
    H_nlos_buffer(ctf_am_reflection, :) = Total_nH_los;
end


%% 1.1.2 
fc_buffer1 = zeros(length(all_reff_num), 2); 
figure;
t_2 = tiledlayout(2,length(all_reff_num));
% title(t_2, "Channel Freq. Corr. and empirical coherence bandwidth")
title(t_2, " Channel Freq. Correlation & Coherence Bandwidth")
% % Enable grid lines


for ctf_am_reflection = 1:length(all_reff_num)

    %For LOS 
    H_current_general = H_los_buffer(ctf_am_reflection,:);
    R_los = zeros(1,floor(n_datapoint/2));
    for change_freq = 1:floor(n_datapoint/2)
        index = 1;
        data_quantity = 0;
        Freq_corr_integr = 0;
        while index + change_freq < n_datapoint
            corre_each = H_current_general(index)*conj(H_current_general(index+change_freq));
            data_quantity = data_quantity + 1;
            Freq_corr_integr = Freq_corr_integr + corre_each;
            index = index + 1;
        end
        R_los(change_freq) = Freq_corr_integr / data_quantity;
    end
    R_los = R_los./ max(abs(R_los)); %normalisation
    nexttile(ctf_am_reflection)
    hold on
   
    title("Frequency Correlation LOS " + all_reff_num(ctf_am_reflection) + " reflection(s)")
    xlabel("\Deltaf [Hz]")
    ylabel("R(\Deltaf)")
plot(w_freq(1:floor(n_datapoint/2)), abs(R_los), 'DisplayName', strcat('los'), 'Color', 'green')
plot(w_freq(1:floor(n_datapoint/2)), horizontal_threshold_fixed, 'DisplayName', strcat('hard threshold'), 'Color', 'magenta')
point = get_threshold(R_los, hard_threshold);
y_min = min(abs(R_los));  % Minimum y-value of the graph
y_max = max(abs(R_los));  % Maximum y-value of the graph

plot((point-1)*w_freq(2)*[1 1], [y_min y_max], ':r', 'DisplayName', strcat('\Deltaf_c =', num2str(w_freq(point)/10^6,3), 'MHz'), 'Color', 'black')

% plot((point-1)*w_freq(2), abs(R_los(point)), 'xr', 'DisplayName', strcat('\Deltaf_c =', num2str(w_freq(point)/10^6,3), 'MHz'), 'Color', 'black')

    legend('show')
    hold off

    fc_buffer1(ctf_am_reflection,1) = point*w_freq(2);

    %For  NLOS
    H_current_general = H_nlos_buffer(ctf_am_reflection,:);
    R_nlos = zeros(1,floor(n_datapoint/2));
    for change_freq = 1:floor(n_datapoint/2)
        index = 1;
        data_quantity = 0;
        Freq_corr_integr = 0;
        while index + change_freq < n_datapoint
            corre_each = H_current_general(index)*conj(H_current_general(index+change_freq));
            data_quantity = data_quantity + 1;
            Freq_corr_integr = Freq_corr_integr + corre_each;
            index = index + 1;
        end
        R_nlos(change_freq) = Freq_corr_integr / data_quantity;
    end
    R_nlos = R_nlos./ max(abs(R_nlos)); % Normalise  the R_nlos for  a better  result


    nexttile(length(all_reff_num) + ctf_am_reflection)
    hold on
   
    title("Frequency Correlation NLOS " + all_reff_num(ctf_am_reflection) + " reflection(s)")
    xlabel("\Deltaf [Hz]")
    ylabel("R(\Deltaf)")
    plot(w_freq(1:floor(n_datapoint/2)),abs(R_nlos),'DisplayName',strcat('nlos'), 'Color', 'green')
    plot(w_freq(1:floor(n_datapoint/2)),horizontal_threshold_fixed,'DisplayName',strcat('hard threshold'),'Color', 'magenta')

    point = get_threshold(R_nlos, hard_threshold);
    % plot(w_freq(2)*(point-1), abs(R_nlos(point)), 'xr', 'DisplayName',strcat('\Deltaf_c =', num2str(w_freq(point)/10^6,3), 'MHz'), 'Color', 'black')

    y_min = min(abs(R_nlos));  % Minimum y-value of the graph
    y_max = max(abs(R_nlos));  % Maximum y-value of the graph

    plot((point-1)*w_freq(2)*[1 1], [y_min y_max], ':r', 'DisplayName', strcat('\Deltaf_c =', num2str(w_freq(point)/10^6,3), 'MHz'), 'Color', 'black')
    legend('show')
    hold off

    fc_buffer1(ctf_am_reflection,2) = point*w_freq(2);
end

%% 1.1.3 Channel Impulse Response
% Deduce the Channel Impulse Response (CIR) and compute the corresponding 
% Power Delay Profile (PDP). Provide reasoning on the obtained plot and the simulated
% multipath channel. Compute the delay spread based on the PDP. Obtain the coherence
% bandwidth through the delay spread and compare it with the one obtained in the previous step.
figure;
tile_4= tiledlayout(2,length(all_reff_num));
title(tile_4,"CIR LOS and NLOS Reflection")

figure;
tiled_5=tiledlayout(2,length(all_reff_num));
title (tiled_5, "PDP LOS AND  NLOS")  %,  " Peak to Hardcut=" + num2str (ch.path_params.max_delay) +  " sec ");

% title (tiled_5, "PDP LOS",  " Peak to Hardcut=" + num2str (ch.path_params.max_delay) +  " sec ");


LOS_peak_store = zeros(length(all_reff_num), 1);
sigma_tau_store = zeros(length(all_reff_num), 2);  %Delay spread value
fc_store_2 = zeros(length(all_reff_num), 2);
P_store = zeros(length(all_reff_num), 2);   



% CIR For LOS and  NLOS
for i= 1: length(all_reff_num)

   max_bin=150;
   max_time_bin=1/ch.B;
   index4_max_delay = floor(ch.path_params.max_delay/max_time_bin);
   tim=(0:n_datapoint-1)*max_time_bin;

   ctf_am_reflection = all_reff_num(i);

   H_current_general = H_los_buffer(i, :);
   H_current_general_time = ifft(H_current_general);
  %  the  plot  of  Channel Impulse Response
    nexttile(tile_4, i)
    title("LOS " +  ctf_am_reflection  + " Reflection(s)")

    hold on
    t_tem=max_bin*max_time_bin;
    % Find the indices where tim is within the desired time range
    valid_indices = tim <= t_tem;

    % Plot only the data within the desired time range
    stem(tim(valid_indices), 10 * log10(abs(H_current_general_time(valid_indices))), 'BaseValue', -80, 'Color', 'k')
    xlabel("Time (s)")
    ylabel("10log_{10}(CIR) (dB)")


    %% Power Delay Profile for LOS
    % calculating the squared magnitude of each sample in the impulse response.
    Pow_Delay_Profil=Power_Delay_Profile(H_current_general_time);
    [~,  LOS_count_position] = max(Pow_Delay_Profil);
    Pow_Delay_Profil = Pow_Delay_Profil(LOS_count_position: index4_max_delay);
    t_temp=tim(LOS_count_position : index4_max_delay);

    nexttile(tiled_5, i)
    title("LOS," +  ctf_am_reflection + " Reflection(s)")

    hold on
    stem(t_temp, 10*log10(Pow_Delay_Profil), 'Basevalue', -150, 'color','b')  % channel attenuation considered
    xlabel("t(s)")
    ylabel("10log_{10}(PDP) (dB)")

    xlim([t_temp(1) t_temp(end)])       % limit  the  signal  along  x axis
    ylim([-150 -90])
 
    % Use  the formula for PDP
    P = Total_Power(Pow_Delay_Profil , max_time_bin );
    Tau_m = Mean_Delay (Pow_Delay_Profil, P, max_time_bin);
    sigma_tau = calculate_delay_spread(Pow_Delay_Profil,  P, Tau_m, max_time_bin); % Delay spread  for LOS
    
    delta_f_c = 1/(2* 3.142* sigma_tau);
    sigma_tau_store (i,1) = delta_f_c;   %storage  of  delay spread value for  table formation
    fc_store_2 (i,2)= delta_f_c;
    P_store (i, 1)= P;


    % Non_Line of Sight
   H_current_general = H_nlos_buffer(i, :);
   H_current_general_time = ifft(H_current_general);

   %  the  plot  of  Channel Impulse Response
    nexttile(tile_4, length(all_reff_num)  + i)
    title("NLOS, " +  ctf_am_reflection + " Reflection(s)")

    
    hold on

    t_tem = max_bin*max_time_bin;
    % Find the indices where tim is within the desired time range
    valid_indices = tim <= t_tem;

    stem(tim(valid_indices), 10*log10(abs(H_current_general_time(valid_indices))), 'BaseValue', -80, 'color','k')

    xlabel("time(s)")
    ylabel("10log_{10}(CIR) (dB)")


    % Power Delay Profile for LOS
    Pow_Delay_Profil=Power_Delay_Profile(H_current_general_time);

    % [~,  LOS_count_position] = max(Pow_Delay_Profil);
    Pow_Delay_Profil = Pow_Delay_Profil(LOS_count_position: index4_max_delay);
    t_temp=tim(LOS_count_position : index4_max_delay);

    nexttile(tiled_5, length(all_reff_num) + i)
    title("NLOS, " +  ctf_am_reflection + " Reflection(s)")

    hold on
    stem(t_temp, 10*log10(Pow_Delay_Profil), 'Basevalue', -150,'color','b')  % channel attenuation considered
    xlabel("t(s)")
    ylabel("10log_{10}(PDP) (dB)")

    xlim([t_temp(1) t_temp(end)])  % limit  the  signal  along  x axis
    ylim([-150 -90])
 
    % Use  the formula for PDP
    P = Total_Power(Pow_Delay_Profil , max_time_bin );
    Tau_m = Mean_Delay (Pow_Delay_Profil, P, max_time_bin);
    sigma_tau = calculate_delay_spread(Pow_Delay_Profil,  P, Tau_m, max_time_bin);
    
    delta_f_c = 1/(2* 3.142* sigma_tau);
    sigma_tau_store (i,2) = sigma_tau;   % store  delay spread  value  for filling NLOS  table
    sigma_tau_store (i,2) = sigma_tau;   % store  delay spread  value  for filling NLOS  table
    fc_store_2 (i,2)= delta_f_c;
    P_store (i, 2)= P;

end

% funcs
function point = get_threshold(array, hard_threshold)
    for avg_ind = 1:length(array)
        if (abs(array(avg_ind)) < hard_threshold)
            point = avg_ind-1;
            return
        end
    end
    point = length(array);
    return
end

%% 
%  Other functions
function PDP = Power_Delay_Profile(h_impulse)
    PDP = zeros(1,length(h_impulse));
    for bin = 1:length(h_impulse)
        PDP(bin) = abs(h_impulse(bin)^2);
    end
end

function P = Total_Power(PDP, Delaytau)
    Delaytau_axis = Delaytau * (0:length(PDP)-1);
    P = trapz(Delaytau_axis, PDP);
end

%The Mean Delay (Tau_m) is the first moment of the PDP, which gives an
%estimate of the average time delay of the channel.
% Mean delay is calculated as the weighted average of the delay bins, where the weights are 
% the power values in each bin. It represents the average time delay experienced by the received signal.
function tau_m = Mean_Delay(PDP, Total_Power, Delaytau)
    Delaytau_axis = Delaytau * (0:length(PDP)-1);
    max_delay_tau = Delaytau_axis.*PDP;
    tau_m = trapz(Delaytau_axis, max_delay_tau) / Total_Power;
end

function Sigma_Tau = calculate_delay_spread(PDP, P, tau_m, Delaytau)
% The Delay Spread (Sigma_Tau) is a measure of the spread in time delays
% experienced by the channel. It is the square root of the second central moment of the PDP, which gives an estimate of the temporal spread of the channel.
% and it quantifies the time spread of the channel's multipath components.
    delay_axis = Delaytau * (0:length(PDP)-1);
    max_delay_tau = (delay_axis.^2).*PDP;
    Sigma_Tau = abs(trapz(delay_axis, max_delay_tau))/P;
    Sigma_Tau = Sigma_Tau - tau_m^2;
    Sigma_Tau = sqrt(Sigma_Tau);
end






