clc; clear;
close all;
addpath functions
%% Appendix parameters
% base station at origin (0.0)
ch.Ptx = 0.1;                                           % transmit power
ch.wt = 8;                                              % top wall (horizontal plane)                     
ch.wb = -8;                                             % bottom wall (horizontal plane)                     
ch.rho = 3;                                             % max number of reflections                     
ch.fc = 3.6e9;                                          % carrier frequency                     
ch.B = 100e6;                                           % system bandwidth                    
ch.Q = 1024;                                            % number of frequency bins/OFDM subcarriers                   
lambda_fc = 3e8/ch.fc;                                  % carrier wavelength                    
dr = lambda_fc/4;                                       % antenna spacing for MIMO-only (you can change the factor 4)                     
ch.beta = 2*pi/lambda_fc;                               % carrier wavenumber                     
ch.A = ch.Ptx*lambda_fc^2/(4*pi)^2;                     % pathloss constant                     
ch.path_params.Npaths = 8;                              % number of cluster paths per main path                     
ch.path_params.max_delay = 500e-9;                      % max path delay, hard cut                     
ch.Rx_pos_x = 12;                                       % Receiver position (x coordinate)
ch.Rx_pos_y = 6;                                        % Receiver position (y coordinate)

%% OFDM

[H_los, H_nlos, RandTheta] = getWideBand(ch); % simulates wideband siso channel


% OFDM parameters
modulationOrder=16; %  number of symbols) for each subcarrier in the OFDM signal.
numOfSymbols=300;          % OFDM symbols that will be transmitted.                                               
numOfSubcarriers = 1024;                    %                   total number of subcarriers used in the OFDM signal.                               
totalSamples=numOfSymbols*numOfSubcarriers;
cyclicPrefixLength = 1024/8;                                   %  1/8th of the length of the symbol, which means that the cyclic prefix will be 128 samples long.
useBPSKModulation = true;                                      %  flag that determines whether  (BPSK) will be used for each subcarrier.
bitErrorRate=[];
LOSChannelImpulseResponse=ifft(H_los);
[~,ind]= max(LOSChannelImpulseResponse);
LOSChannelImpulseResponse=LOSChannelImpulseResponse(ind:ind+50); %portion of the LOS channel impulse response that includes the main lob

% LOSChannelImpulseResponse=[1];
H_los = fft(LOSChannelImpulseResponse, ch.Q);
NLOSChannelImpulseResponse=ifft(H_nlos);
errorVector=[]; %store the error rate for each CFO value tested.

numOfPreambles=2;
estimateCFO=true; % flag that determines whether carrier frequency offset (CFO) will be estimated and corrected.
CFOVector = [0.01 0.1 1 5 10 15]; %CFO values that will be tested in the simulatio
errorVector=zeros(length(CFOVector),11);
berMat=zeros(length(CFOVector), 11);
%figure;

% Call the function to calculate the bitErrorRate
for i=1:length(CFOVector)
    cfo_ppm=CFOVector(i)
    [bitErrorRate,err, signalToNoiseRatio, errorVector(i,:)] = Freq_ofdm_acq(H_los,true,1024,300, cfo_ppm,ch,50);
    beravg=mean(bitErrorRate);                           % Calculate the average bitErrorRate over multiple trials
    berMat(i,:)=beravg;                                  % Store the average bitErrorRate value
    semilogy((beravg))                                   % Plot the bitErrorRate vs signalToNoiseRatio curve for the current CFO value
    hold on
end
legend("CFO=0.01","CFO=0.1","CFO=1","CFO=5","CFO=10","CFO=15")
xlabel("SNR")
ylabel("BER")
title('Bit Error Rate FOR DIFFERENT CFO')


function [averageBERStore,errorValue, signalToNoiseRatio, errorVector] = Freq_ofdm_acq(H,estimateCFO,numOfSubcarriers,numOfSymbols, cfo_ppm, ch, numOfExperiments)
%% OFDM

am_preamble=2;
SNR_db_collection= 0:2:40;

averageBERStore=zeros(numOfExperiments,length(SNR_db_collection));

modulationOrder=16; % 

totalSamples=numOfSymbols*numOfSubcarriers;

cyclicPrefixLength = 1024/8;    % 1/8th of the length of the symbol, which means that the cyclic prefix will be 128 samples long.

channelImpulseResponse=ifft(H);

[~,ind]= max(channelImpulseResponse);
channelImpulseResponse=channelImpulseResponse(ind:ind+50);
H = fft(channelImpulseResponse, numOfSubcarriers);

for k=1:numOfExperiments
    k
    errorVector=zeros(1,6 ); %length(0:4:40)
    bitErrorRate=zeros(1,6); 
for i=1:length(SNR_db_collection)
signalToNoiseRatio=SNR_db_collection(i);


% Generate random dataBits bits
dataBits = randsrc(1,totalSamples, 0:modulationOrder-1);

% Calculate the average power of the QAM symbols
avgPower = mean(abs(qammod(0:modulationOrder-1, modulationOrder)).^2);

% Normalize the modulated symbols to have unit average power
transmittedData = qammod(dataBits, modulationOrder) / sqrt(avgPower);

                

% dataBits = randsrc(1,totalSamples, 0:modulationOrder-1);
% transmittedData=qammod(dataBits,modulationOrder, "UnitAveragePower", 1);

if estimateCFO
    preambleData = 2*randi([0 1], 1, numOfSubcarriers)-1;

%     preambleData=randsrc(1,numOfPreambles*numOfSubcarriers, 0:modulationOrder-1);
    transmittedData = [preambleData preambleData transmittedData];
    numOfSymbols = numOfSymbols + am_preamble;
    totalSamples = numOfSymbols*numOfSubcarriers;
end


transmittedDataMatrixFrequency=reshape(transmittedData,numOfSubcarriers,numOfSymbols);
transmittedDataMatrixTime=ifft(transmittedDataMatrixFrequency, numOfSubcarriers,1);



% Generate OFDM symbol
transmittedDataMatrixTimeWithPrefix = [transmittedDataMatrixTime(end - cyclicPrefixLength+1:end,:); transmittedDataMatrixTime];
transmittedOFDMData=transmittedDataMatrixTimeWithPrefix(:);

transmittedOFDMData=conv(transmittedOFDMData, channelImpulseResponse);

transmittedOFDMData=transmittedOFDMData(1:numel(transmittedDataMatrixTimeWithPrefix));

% Add noise to OFDM symbol receivedSignal=ofdm_symbol;
timeIndex = 0:length(transmittedOFDMData)-1;
    
       

carrierFrequencyOffset=exp(-1j*(2*pi*ch.fc*cfo_ppm*1e-6)*timeIndex.'./ch.B);
transmittedOFDMData = carrierFrequencyOffset.*transmittedOFDMData;
rx_data_fft_wpf = awgn(transmittedOFDMData, signalToNoiseRatio, 'measured');

% __________________________________________

rx_data_fft_with_prefix=reshape(rx_data_fft_wpf,numOfSubcarriers+cyclicPrefixLength,numOfSymbols);

% Remove cyclic prefix
f_cfo_est = angle(rx_data_fft_with_prefix(:, 1).' * conj(rx_data_fft_with_prefix(:, 2)))/(2*pi*(ch.Q+cyclicPrefixLength)/ch.B);
errorVector(i)=(f_cfo_est-ch.fc*cfo_ppm*1e-6);

%%  CFO Correction
correction=exp(1j*2*pi*f_cfo_est*timeIndex.'./ch.B);
rx_data_fft_with_prefix=correction.*rx_data_fft_wpf;
rx_data_fft_with_prefix=reshape(rx_data_fft_with_prefix,numOfSubcarriers+cyclicPrefixLength,numOfSymbols);


% Perform FFT on received dataBits
fft_data = fft(rx_data_fft_with_prefix(cyclicPrefixLength+1:end, :), numOfSubcarriers,1);

if estimateCFO  
    rx_serial_fft_data=reshape(fft_data,1,totalSamples);

    preamble_rx = rx_serial_fft_data(1:numOfSubcarriers);

    preamble_original=transmittedData(1:numOfSubcarriers);

    H_est=preamble_rx./preamble_original;

    preamble_rx = rx_serial_fft_data(numOfSubcarriers+1:2*numOfSubcarriers);

    preamble_original=transmittedData(numOfSubcarriers+1:2*numOfSubcarriers);

    H_est=H_est+(preamble_rx./preamble_original);
    H_est=H_est/2;

else
    H_est=H;
end
    
if estimateCFO
    fft_data=fft_data(:);
     for v = 1:numOfSymbols
               rx_data_equalized((v-1)*numOfSubcarriers+1:v*numOfSubcarriers) = fft_data((v-1)*numOfSubcarriers+1:v*numOfSubcarriers)./H_est.'; %ZF equalizer
     end
else
    rx_data_equalized = fft_data ./ H.';
end
errorValue=mse(H,H_est);

rx_data_demod = qamdemod(rx_data_equalized(2*numOfSubcarriers+1:end),modulationOrder, "UnitAveragePower", 1 );

% Calculate bit errorVector rate

bitErrorRate(i) = (sum(de2bi(rx_data_demod) ~= de2bi(dataBits), "all")) / numel(de2bi(dataBits));

end
averageBERStore(k,:)=bitErrorRate;
end
end





















