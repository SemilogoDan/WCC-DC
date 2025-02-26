addpath functions
clear all
%%
ch.Ptx = 0.1;                                 % transmit power
ch.wt = 8;                                    % top wall
ch.wb = -8;                                   % bottom wall
ch.rho = 3;                                   % max number of reflections
ch.fc = 3.6e9;                                % carrier frequency
ch.B = 100e6;                                 % system bandwidth
ch.Q = 1024;                                  % number of frequency bins/OFDM subcarriers
lambda_fc = 3e8/ch.fc;                        % carrier wavelength
dr = lambda_fc/4;                             % antenna spa cing for MIMO-only (you can change the factor 4)
ch.beta = 2*pi/lambda_fc;                     % carrier wavenumber
ch.A = ch.Ptx*lambda_fc^2/(4*pi)^2;           % pathloss constant
ch.path_params.Npaths = 8;                    % number of cluster paths per main path
ch.path_params.max_delay = 500e-9;            % max path delay, hard cut

%% WideBand
ch.Rx_pos_x = 12;
ch.Rx_pos_y = 6;

[H_los, H_nlos, RandTheta] = getWideBand(ch); % simulates wideband siso channel

% Call the modified OFDM function with channel estimation and BER evaluation
[ber1, LOS_err1, snr1] = ofdm(H_los, false, ch.Q, 40); % With channel estimation
[ber3, NLOS_err3, snr2] = ofdm(H_nlos, false, ch.Q, 40); % With channel estimation


% Create a customized BER vs. SNR plot
figure;
semilogy(-4:4:40, ber1, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'LOS');
hold on;
semilogy(-4:4:40, ber3, '-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'NLOS');
title('BER OFDM Transceiver without Channel Estimation');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
grid on;
legend('Location', 'best');
legend('show');
% 
% disp(['Channel estimation Error LOS: ',num2str(LOS_err1)])
% disp(['Channel estimation Error NLOS: ',num2str(NLOS_err3)])
% disp(['Channel estimation Error LOS: ', num2str(abs(LOS_err1))])
% disp(['Channel estimation Error NLOS: ', num2str(abs(NLOS_err3))])
% Display channel estimation errors
disp(['Channel estimation Error LOS: ', num2str(LOS_err1, '%.15f')]);
disp(['Channel estimation Error NLOS: ', num2str(NLOS_err3, '%.15f')]);
if LOS_err1>NLOS_err3
    disp('LOS worse has ber than NLOS')
else
    disp('NLOS worse has ber than LOS')
end

function [ber,MSE_err, snr] = ofdm(H,estimate,subcarrier_amt,symbol_amt)
%% OFDM
close all

%____________
% OFDM parameters
M=16; % 
tot_samples=symbol_amt*subcarrier_amt;
cp_len = 64;    % Cyclic prefix length
ber=[];
CIR=ifft(H);

[~,ind]= max(CIR);
CIR=CIR(ind:ind+50);
H = fft(CIR, subcarrier_amt);


preamble_amt=2;


for snr= -4:4:40

% Generate random data bits
data = randsrc(1,tot_samples, 0:M-1);
data_tx=qammod(data,M, "UnitAveragePower", 1);

if estimate
    preamble = 2*randi([0 1], 1, subcarrier_amt)-1;

%     preamble=randsrc(1,preamble_amt*subcarrier_amt, 0:M-1);
    data_tx = [preamble preamble data_tx];
    symbol_amt = symbol_amt + preamble_amt;
    tot_samples = symbol_amt*subcarrier_amt;
end


data_tx_mat_freq=reshape(data_tx,subcarrier_amt,symbol_amt);
data_tx_mat_time=ifft(data_tx_mat_freq, subcarrier_amt,1);



% Generate OFDM symbol
data_tx_mat_time_wpf = [data_tx_mat_time(end-cp_len+1:end,:); data_tx_mat_time];
ofdm_data_tx=data_tx_mat_time_wpf(:);

ofdm_data_tx=conv(ofdm_data_tx, CIR);

ofdm_data_tx=ofdm_data_tx(1:numel(data_tx_mat_time_wpf));

% Add noise to OFDM symbol
rx_data_fft_wpf = awgn(ofdm_data_tx, snr, 'measured');

% __________________________________________

rx_data_fft_wpf=reshape(rx_data_fft_wpf,subcarrier_amt+cp_len,symbol_amt);
% Remove cyclic prefix
rx_data_fft = rx_data_fft_wpf(cp_len+1:end,:);

% Perform FFT on received data
% subcarrier_amt=
subcarrier_amt = 1024;         % Number of subcarriers
fft_data = fft(rx_data_fft, subcarrier_amt,1);
% rx_serial_fft_data=reshape(fft_data,1,tot_samples);

if estimate  
    rx_serial_fft_data=reshape(fft_data,1,tot_samples);
    preamble_rx = rx_serial_fft_data(1:subcarrier_amt);
    preamble_original=data_tx(1:subcarrier_amt);
    H_est=preamble_rx./preamble_original;

    preamble_rx = rx_serial_fft_data(subcarrier_amt+1:2*subcarrier_amt);
    preamble_original=data_tx(subcarrier_amt+1:2*subcarrier_amt);
    H_est=H_est+(preamble_rx./preamble_original);
    
    H_est=H_est/2;
else
    H_est=H;
end


if estimate
    rx_data_equalized = fft_data(:,3:end) ./ H.';
else
    fft_data=fft_data(:);  %ZF equalization 
     for v = 1:symbol_amt
               rx_data_equalized((v-1)*subcarrier_amt+1:v*subcarrier_amt) = fft_data((v-1)*subcarrier_amt+1:v*subcarrier_amt)./H_est.'; %ZF equalizer
     end
end
MSE_err=mse(H,H_est)

% Demodulate received data
rx_data_demod = qamdemod(rx_data_equalized(:),M, "UnitAveragePower", 1 );

% Calculate bit error rate

ber = [ber (sum(de2bi(rx_data_demod) ~= de2bi(data), "all")) / numel(de2bi(data))];


% fprintf('Bit error rate: %f\n', ber);
end

end

% addpath functions
% clear all
% %%
% ch.Ptx = 0.1;                                 % transmit power
% ch.wt = 8;                                    % top wall
% ch.wb = -8;                                   % bottom wall
% ch.rho = 3;                                   % max number of reflections
% ch.fc = 3.6e9;                                % carrier frequency
% ch.B = 100e6;                                 % system bandwidth
% ch.Q = 1024;                                  % number of frequency bins/OFDM subcarriers
% lambda_fc = 3e8/ch.fc;                        % carrier wavelength
% dr = lambda_fc/4;                             % antenna spa cing for MIMO-only (you can change the factor 4)
% ch.beta = 2*pi/lambda_fc;                     % carrier wavenumber
% ch.A = ch.Ptx*lambda_fc^2/(4*pi)^2;           % pathloss constant
% ch.path_params.Npaths = 8;                    % number of cluster paths per main path
% ch.path_params.max_delay = 500e-9;            % max path delay, hard cut
% 
% %% WideBand
% ch.Rx_pos_x = 12;
% ch.Rx_pos_y = 6;
% 
% [H_los, H_nlos, RandTheta] = getWideBand(ch); % simulates wideband siso channel
% 
% % Call the modified OFDM function with channel estimation and BER evaluation
% [ber1, LOS_err1, snr1] = ofdm(H_los, false, ch.Q, 40); % With channel estimation
% [ber3, NLOS_err3, snr2] = ofdm(H_nlos, false, ch.Q, 40); % With channel estimation
% 
% 
% % Create a customized BER vs. SNR plot
% figure;
% semilogy(-4:4:40, ber1, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'LOS');
% hold on;
% semilogy(-4:4:40, ber3, '-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'NLOS');
% title('BER OFDM Transceiver without Channel Estimation');
% xlabel('SNR (dB)');
% ylabel('Bit Error Rate (BER)');
% grid on;
% legend('Location', 'best');
% legend('show');
% % 
% % disp(['Channel estimation Error LOS: ',num2str(LOS_err1)])
% % disp(['Channel estimation Error NLOS: ',num2str(NLOS_err3)])
% % disp(['Channel estimation Error LOS: ', num2str(abs(LOS_err1))])
% % disp(['Channel estimation Error NLOS: ', num2str(abs(NLOS_err3))])
% % Display channel estimation errors
% disp(['Channel estimation Error LOS: ', num2str(LOS_err1, '%.15f')]);
% disp(['Channel estimation Error NLOS: ', num2str(NLOS_err3, '%.15f')]);
% if LOS_err1>NLOS_err3
%     disp('LOS worse has ber than NLOS')
% else
%     disp('NLOS worse has ber than LOS')
% end
% 
% function [ber,MSE_err, snr] = ofdm(H,estimate,subcarrier_amt,symbol_amt)
% %% OFDM
% close all
% 
% %____________
% % OFDM parameters
% M=16; % 
% tot_samples=symbol_amt*subcarrier_amt;
% cp_len = 64;    % Cyclic prefix length
% ber=[];
% CIR=ifft(H);
% [~,ind]= max(CIR);
% CIR=CIR(ind:ind+50);
% H = fft(CIR, subcarrier_amt);
% 
% 
% preamble_amt=2;
% 
% 
% for snr= -4:4:40
% 
% % Generate random data bits
% data = randsrc(1,tot_samples, 0:M-1);
% data_tx=qammod(data,M, "UnitAveragePower", 1);
% 
% if estimate
%     preamble = 2*randi([0 1], 1, subcarrier_amt)-1;
% 
% %     preamble=randsrc(1,preamble_amt*subcarrier_amt, 0:M-1);
%     data_tx = [preamble preamble data_tx];
%     symbol_amt = symbol_amt + preamble_amt;
%     tot_samples = symbol_amt*subcarrier_amt;
% end
% 
% 
% data_tx_mat_freq=reshape(data_tx,subcarrier_amt,symbol_amt);
% data_tx_mat_time=ifft(data_tx_mat_freq, subcarrier_amt,1);
% 
% 
% 
% % Generate OFDM symbol
% data_tx_mat_time_wpf = [data_tx_mat_time(end-cp_len+1:end,:); data_tx_mat_time];
% ofdm_data_tx=data_tx_mat_time_wpf(:);
% 
% ofdm_data_tx=conv(ofdm_data_tx, CIR);
% 
% ofdm_data_tx=ofdm_data_tx(1:numel(data_tx_mat_time_wpf));
% 
% % Add noise to OFDM symbol
% rx_data_fft_wpf = awgn(ofdm_data_tx, snr, 'measured');
% 
% % __________________________________________
% 
% rx_data_fft_wpf=reshape(rx_data_fft_wpf,subcarrier_amt+cp_len,symbol_amt);
% % Remove cyclic prefix
% rx_data_fft = rx_data_fft_wpf(cp_len+1:end,:);
% 
% % Perform FFT on received data
% % subcarrier_amt=
% subcarrier_amt = 1024;         % Number of subcarriers
% fft_data = fft(rx_data_fft, subcarrier_amt,1);
% % rx_serial_fft_data=reshape(fft_data,1,tot_samples);
% 
% if estimate  
%     rx_serial_fft_data=reshape(fft_data,1,tot_samples);
%     preamble_rx = rx_serial_fft_data(1:subcarrier_amt);
%     preamble_original=data_tx(1:subcarrier_amt);
%     H_est=preamble_rx./preamble_original;
%     preamble_rx = rx_serial_fft_data(subcarrier_amt+1:2*subcarrier_amt);
%     preamble_original=data_tx(subcarrier_amt+1:2*subcarrier_amt);
%     H_est=H_est+(preamble_rx./preamble_original);
%     H_est=H_est/2;
% else
%     H_est=H;
% end
% 
% 
% if estimate
%     rx_data_equalized = fft_data(:,3:end) ./ H.';
% else
%     fft_data=fft_data(:);  %ZF equalization 
%      for v = 1:symbol_amt
%                rx_data_equalized((v-1)*subcarrier_amt+1:v*subcarrier_amt) = fft_data((v-1)*subcarrier_amt+1:v*subcarrier_amt)./H_est.'; %ZF equalizer
%      end
% end
% MSE_err=mse(H,H_est)
% 
% % Demodulate received data
% rx_data_demod = qamdemod(rx_data_equalized(:),M, "UnitAveragePower", 1 );
% 
% % Calculate bit error rate
% 
% ber = [ber (sum(de2bi(rx_data_demod) ~= de2bi(data), "all")) / numel(de2bi(data))];
% 
% 
% % fprintf('Bit error rate: %f\n', ber);
% end
% 
% end
% 
