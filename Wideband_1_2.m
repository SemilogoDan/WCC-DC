%% Basic  OFDM  Signal Generation
clc; clear;
close all;
addpath functions
%% Appendix parameters
% base station at origin (0.0)
ch.Ptx = 0.1;                       % transmit power
ch.wt = 8;                          % top wall (horizontal plane)
ch.wb = -8;                         % bottom wall (horizontal plane)
ch.rho = 3;                         % max number of reflections
ch.fc = 3.6e9;                      % carrier frequency
ch.B = 100e6;                       % system bandwidth
ch.Q = 1024;                        % number of frequency bins/OFDM subcarriers
lambda_fc = 3e8/ch.fc;              % carrier wavelength
dr = lambda_fc/4;                   % antenna spacing for MIMO-only (you can change the factor 4)
ch.beta = 2*pi/lambda_fc;           % carrier wavenumber
ch.A = ch.Ptx*lambda_fc^2/(4*pi)^2; % pathloss constant
ch.path_params.Npaths = 8;          % number of cluster paths per main path
ch.path_params.max_delay = 500e-9;  % max path delay, hard cut=5e7
ch.Rx_pos_x = 12;
ch.Rx_pos_y = 6;

cp_length_est = ch.path_params.max_delay*ch.B;   %50
%% 1.2.1 OFDM simulation
% EDIT SIMULATION HERE
am_symbols = 2^16;
bitspersymbol_collection = [2,4,8,16]; %check different values of QAM
SNR_db_collection =-5:2:40;
%SNR_db_collection =-4:0.5:12.5;

LOS = 1;
% cp_length = 0;
cp_length = 2^8; % in symbols per block; set to 0 to disable
am_realizations = 30;
% 2*ch.B represents the total bandwidth allocated for OFDM.
% ch.Q is the number of subcarriers in the OFDM system.
% (2*ch.B)/ch.Q calculates the spacing between adjacent subcarriers.
% (2*ch.B-2*ch.B/ch.Q) sets the last [subcarrier's center frequency' ...
%     ' to be just below the maximum frequency of the allocated bandwidth.]
% %other constants
f = 0:(2*ch.B)/ch.Q:(2*ch.B-2*ch.B/ch.Q); % Obey Nyquist rate, frequency values for each subcarrier's center frequency
n = length(f); %amount of datapoints in sample sequence
time_bin = 1/(2*ch.B); %in [s]: time resolution (time  duration  of OFDM of ofdm symbol)
max_delay_index = floor(ch.path_params.max_delay / time_bin); % 500 time-domain samples for the channel response to account for the maximum delay of 0.5 nanoseconds
t = (0:n-1)*time_bin; % time vector with n element scaled over time_bin



%% OBTAIN BER

BER_store = zeros(length(bitspersymbol_collection),length(SNR_db_collection),am_realizations);

for SNR_db_index = 1:length(SNR_db_collection)
    SNR_db = SNR_db_collection(SNR_db_index);

    for bitspersymbol_index = 1:length(bitspersymbol_collection)
        bitspersymbol = bitspersymbol_collection(bitspersymbol_index);
        for realization_index = 1:am_realizations
              %frequency-domain channel frequency responses for both LOS and NLOS conditions
            % [H_los_standard, H_nlos_standard, RandTheta_standard] = getNarrowBand(ch); %FOR TESTING ONLY%
            [H_los_standard, H_nlos_standard, RandTheta_standard] = getWideBand(ch); % simulates wideband siso channel
            % H_los_standard, H_nlos_standard, getWideBand, ch are all in frequency domain
            h_los_standard = ifft(H_los_standard);

            h_nlos_standard = ifft(H_nlos_standard); %channel's behavior in the time domain,
            % including the effects of multipath propagation and delay spread.

            if LOS == 1
                H_current = H_los_standard;   % freq domain LOS
                h_current = h_los_standard;   % time domain LOS
            else
                H_current = H_nlos_standard;  % freq domain NLOS
                h_current = h_nlos_standard;  %  time domain  NLOS
            end
            %  if LOS == 1
            %     H_current = H_los_standard;   % freq domain LOS
            %     h_current = h_los_standard;   % time domain LOS
            % 
            %     [~, LOS_index] = max(abs(h_current)); % Find the index of the LOS peak
            %     max_delay_index = min(LOS_index + floor(ch.path_params.max_delay / time_bin), length(h_current));
            %     % Ensure max_delay_index does not exceed the length of h_current
            % 
            %     h_current = h_current(LOS_index:max_delay_index); % Extract the channel response within the desired range
            % else
            %     H_current = H_nlos_standard;  % freq domain NLOS
            %     h_current = h_nlos_standard;  %  time domain  NLOS
            % 
            %     % Note: You may need to define the values of LOS_index and max_delay_index for NLOS condition
            %     % They could be calculated in a similar way as for LOS, based on your channel model.
            % end


            %cut from LOS peak until hard cutoff
            if LOS == 1
                [~, LOS_index] = max(h_current);
                h_current = h_current(LOS_index: max_delay_index); %max_delay_index is a given value for hardcut
                H_current = fft(h_current, ch.Q);  % update the  freq domain LOS using the time domain h_current
                %H_current contain complex values, amplitude and phase of
                %the corresponding frequency component.
                % function call that calculates the FFT of the signal h_current using ch.Q points and returns the frequency-domain
            end

            %block lengths
            % am_symbols = 2^16;
            am_bits = am_symbols * bitspersymbol;               % Total number of bits to be transmitted at the transmitter (bits)
            am_blocks = am_symbols/ch.Q;                        %Q the number of frequency bin whihch means the number of subcarriers
            if cp_length ~= 0
                am_bits_padded = am_bits + cp_length*am_blocks*bitspersymbol;
                am_symbols_padded = am_symbols + cp_length*am_blocks;
            else
                am_bits_padded = am_bits;
                am_symbols_padded = am_symbols;
            end

            %get random bit sequence or generate random data bits
            bits = randi([0,1],am_bits,1);

            Symbols = QAM_modulation(bits, bitspersymbol);
    %         qammod()
            symbols = zeros(1,am_symbols);


            if mod(am_symbols, ch.Q) ~= 0
                disp("error: incomplete block-handling has not yet been implemented")
            else
                disp("Data transmission in progress!")
            end

            for i = 1:am_blocks
                begin_index = (i-1)*ch.Q + 1;
                end_index = i*ch.Q;

                % OFDM: IFFT
                symbols(begin_index:end_index) = ifft(Symbols(begin_index:end_index));
            end

            %% add padding by ading  Cyclo Prfix
            % The cycloprefix length is the number of samples copied from the end of an OFDM
            % symbol and inserted at the beginning of the symbol to form the cyclic prefix.
            % It helps in mitigating the effects of inter-symbol interference (ISI) caused by multipath propagation.
            symbols_padded = zeros(1, am_symbols_padded);
            if cp_length == 0
                symbols_padded = symbols;
            else
                for i = 1:am_blocks
                     begin_index = (i-1)*ch.Q + 1;
                     end_index = i*ch.Q;   
                     symbols_padded_start = (i-1)*(ch.Q + cp_length) + 1;
                     symbols_padded_end = (i)*(ch.Q + cp_length);
                     symbols_padded(symbols_padded_start:symbols_padded_start + cp_length - 1) = symbols(end_index-cp_length + 1: end_index);
                     symbols_padded(symbols_padded_start + cp_length : symbols_padded_end) = symbols(begin_index: end_index);
                end
            end

            % OFDM: TOTAL convolution with the channel
            % ofdm_data_tx=conv(ofdm_data_tx, CIR);
            %symbols_padded= symbol transmitted
            % h_current= CIR
            r_noiseless_padded = conv(symbols_padded, h_current);

            %scale it (ONLY FOR PLOT)
            maxvalue_sym = max(abs(symbols_padded));
            maxindex = find(abs(symbols) == maxvalue_sym);
            maxindex = maxindex(1);
            maxvalue_r_noiseless = abs(r_noiseless_padded(maxindex));
            r_noiseless_scaled_padded = r_noiseless_padded*maxvalue_sym/maxvalue_r_noiseless;
            
            % add noise
            % r_noiseless_padded= OFDM data transmitted
            received_padded_signa = awgn(r_noiseless_padded,SNR_db, 'measured');
            %r = r_noiseless; %for testing only!
            r_scaled_padded = received_padded_signa*maxvalue_sym/maxvalue_r_noiseless;



            % cut off cyclic prefix + tail  
            if(length(received_padded_signa) ~= length(symbols_padded))
                boot_cutoff_length = length(h_current);   %h_current= channel impulse response 
                r_padded_no_boot = received_padded_signa(1:end-boot_cutoff_length+1); %removes the cyclic prefix from the received signal discarding the initial boot_cutoff_length samples
            else
                r_padded_no_boot = received_padded_signa; % if no cyclic prefix in the received signa
            end
            Remov_CP_r_symbol_ = zeros(size(symbols));
            for i = 1:am_blocks
                begin_index = (i-1)*ch.Q + 1;
                end_index = i*ch.Q;   
                symbols_padded_start = (i-1)*(ch.Q + cp_length) + 1;
                symbols_padded_end = (i)*(ch.Q + cp_length);
                % copied  to r array excluding the cyclic prefix
                Remov_CP_r_symbol_(begin_index:end_index) = r_padded_no_boot(symbols_padded_start + cp_length:symbols_padded_end);
            end
           % Each column of this matrix represents one block of received OFDM symbols, including the cyclic prefix
            r_temp = reshape(r_padded_no_boot, ch.Q+cp_length, am_blocks);
 
           %Remove  cyclo prefix
            Remov_CP_r_symbol_ = r_temp(cp_length+1:end, :); %effectively removes the cyclic prefix from each block of received symbols..
            %Perform FFT on  received data

            rx_data_equalized = fft(Remov_CP_r_symbol_, ch.Q, 1);
            % Perform  channel  equalization for distorted freq
            % compensation caused  by the channel
            I_est = rx_data_equalized./H_current.';
            I_est = I_est(:); %  convert  from 2 dimensional vector  to 1 dim
            % err= mse(H_current, )
            % Perform channel equalization for distorted frequency compensation caused by the channel
    
           % Calculate Mean Squared Error (MSE)


           % bitspersymbol=M
           %Demodulate  received  data
            bits_restored = QAM_demodulation(I_est,bitspersymbol);
            BER = 0;
            for i = 1:am_bits
                if bits(i) ~= bits_restored(i)
                    BER = BER + 1;
                end
            end
            % calculate bit error rate
            BER = BER / am_bits; %high!
            BER_store(bitspersymbol_index,SNR_db_index, realization_index) = BER;
        end
    end
end
BER_store = mean(BER_store,3);  % MSE
% first dimension corresponds to different 
% modulation schemes, and the second dimension corresponds to different SNR values.
%% plot results
figure;
hold on;
colors = lines(length(bitspersymbol_collection));  % Define colors for each line
markers = {'o', 's', 'd', 'v'};  % Define markers for each line
for bitspersymbol_index = 1:length(bitspersymbol_collection)
    plot(SNR_db_collection, BER_store(bitspersymbol_index, :), ...
        'DisplayName', strcat('Nbps=', num2str(bitspersymbol_collection(bitspersymbol_index))), ...
        'Color', colors(bitspersymbol_index, :), ...
        'Marker', markers{bitspersymbol_index});
end
set(gca, 'YScale', 'log');
grid on;
% title("BER of OFDM with AWGN (CP=0)");
title("BER  VS SNR over Wideband OFDM cp=" + num2str(cp_length))
xlabel("SNR [dB]");
ylabel("BER");
legend('show');
















