clear; close all; clear all;
clc;
%% Initialization
SNR_list    = -10:20:30;       % SNR grid
Nbps        = 2;               % Number of bits per quam symbol
% OFDM parameters
B          = 100e6;            % Bandwidth
T=1/(B);
N_subcrr   = 1024;             % Number of subcarriers
Lcp        = 64;               % CP length, samples (TD)
preamble_L = 2;                % Preamble length, OFDM symbols (FD)
data_L     = 30;               % Data length in # symbols per channel (FD)
f_dc       = 3.6e9;            % Carrier frequency


% some extra controls
% piloting
pilot_enable = 1; %1 = true, 0 = false
pilot_spacing = 128; %we only perform piloting every %pilot_spacing% carriers
%t_0
t_0_enable = 1;

%  QAM Modulation..
Nbits = data_L * N_subcrr * Nbps; %bits
bits = randi([0 1], [1 Nbits]);
bits_tx = bits' ;

%rule we made up
 if(Nbps > 1)
    modulation_mode = 'qam';   
 else
    modulation_mode = 'pam';
 end
  
 [symb_tx] = mapping(bits_tx,Nbps,modulation_mode);

% OFDM Transmitter: 
S2P_tx = reshape(symb_tx,N_subcrr,data_L);%1024 x 300
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Add the preamble
% Nbps=4;
% % BPSK: generate a series of bits and use it for by premble by its BPSK  formation
% bits = randi([0 1], [1 (N_subcrr)*4]);
% bits_tx = bits' ;
% preamble = mapping(bits_tx,Nbps,'qam');
% save('preamble.mat','preamble');
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if(t_0_enable)
    load('preamble.mat');
    preamble_tx=[preamble preamble];
    S2P_tx=[preamble_tx S2P_tx];
end

% Pilot addition
if(pilot_enable)
     pilot_bits = ones(Nbps*2,1);
     pilot_symbol = mapping(pilot_bits,Nbps,modulation_mode); % From MA1 Prof. Francois Modulation and Coding lecture
     if(t_0_enable)
         t_0_offset = 2;
     else
         t_0_offset = 0;
     end
     for i=1+t_0_offset:data_L+t_0_offset %skip the first two since they are preamble
        for j=1:pilot_spacing:N_subcrr
            S2P_tx(j,i)=1;
        end
     end
 end

% IFFT
signal_tx_time = ifft(S2P_tx,[],1);
 
% Addition of CP
signal_tx_time_CP=signal_tx_time(end-(Lcp-1):end,:);
signal_tx_time_total=[signal_tx_time_CP;signal_tx_time];
L_td = size(signal_tx_time_total,1);

% Get the vector: signal_tx & The second propagation path
[m,n]=size(signal_tx_time_total);
mul_tx = m*n;
signal_tx =reshape(signal_tx_time_total,1, mul_tx );
length_constant=length(signal_tx);

chtaps = [1 0.1*randn(1)*exp(-1i*pi*randn(1))];
% chtaps = [1 0 0.4*exp(-1i*pi/8)];
equalization__chtaps= (fft(chtaps,N_subcrr))';

% Multipath propagation
signal_tx=Multipath_Propa(signal_tx,chtaps);
% MRC
num_antennas=1;
for k = 1: length(SNR_list)
    if SNR_list(k)>=0
        SNR_list(k)=SNR_list(k)*num_antennas;
    else
        SNR_list(k)=SNR_list(k)/num_antennas;
    end
end

% Channel propagation:
for k = 1: length(SNR_list)
EbN0 = SNR_list(k);% SNR
energy_of_signal=1/B * trapz(abs(signal_tx).^2);
Eb = energy_of_signal / Nbits ;
N0 = Eb/10.^(EbN0/10);
NoisePower = 2*N0*B;
% OFDM Receiver:
%% 3.0 Add the uncertainty on time-of-arrival(STO) & Add noise
signal_STO=add_STO(signal_tx);
noise = sqrt(NoisePower/2)*(randn(length(signal_STO),1)+...
1i*randn(length(signal_STO),1))';
signal_STO_noise=signal_STO+noise;

%% 3.2_1 Add  CFO   to Received Signal
signal_rx_STO_CFO=add_CFO(signal_STO_noise);
%time_acquisition=auto_corr(signal_rx_STO_CFO)

% Auto-correlation for STO & CFO and compensating
[time_acquisition, CFO_est]=auto_corr(signal_rx_STO_CFO)
CFO__est_list(k)=CFO_est;


% SIMO
taps2 = [1 0.1*randn(1)*exp(-1i*pi*randn(1)) ];
taps3 = [1 0.1*randn(1)*exp(-1i*pi*randn(1)) ];
taps4 = [1 0.1*randn(1)*exp(-1i*pi*randn(1)) ];

% SIMO
[equlization_RX2, time_acquisition_RX2, CFO_est_RX2]=SIMO(signal_tx,taps2,preamble_tx);
[equlization_RX3, time_acquisition_RX3, CFO_est_RX3]=SIMO(signal_tx,taps3,preamble_tx);
[equlization_RX4, time_acquisition_RX4, CFO_est_RX4]=SIMO(signal_tx,taps4,preamble_tx);

% take an average for equlization to use later
time_acquisition_total=[time_acquisition time_acquisition_RX2 time_acquisition_RX3 time_acquisition_RX4];
time_acquisition_mean=round(mean(time_acquisition_total));

% take an average for equlization to use later

CFO_est_total=[CFO_est CFO_est_RX2 CFO_est_RX3 CFO_est_RX4];
CFO_est_mean=mean(CFO_est_total);

%% compensating STO CFO
signal_rx=compensating(signal_rx_STO_CFO,CFO_est,T,time_acquisition,length_constant);

% Serial to parallel & Remove of CP
signal_rx_fx=End_Delete_CP(signal_rx,m,n,Lcp);

% Premble and get the channel impulse response
preamble_rx1 = signal_rx_fx(:,1);
preamble_rx2 = signal_rx_fx(:,2);
preamble_rx = [preamble_rx1 preamble_rx2];
CH_F = mean(preamble_rx./preamble_tx,2);
CH_T = ifft(CH_F);
CH_T = [CH_T(1:64);zeros((N_subcrr-Lcp),1)];
CH_F_t = fft(CH_T);

%%  CH_F(Channel estimation in Freq domain)
equalization = CH_F;
signal_rx_est_freq=((signal_rx_fx)./(equalization));
signal_rx_fx=signal_rx_est_freq;
% CH_T
equalization2= CH_F_t;
CIR=ifft(equalization2);% CIR(1:20)
figure;plot(abs(CIR));

%% 3.0_3   MSE calculate
MSE_F(k)=Mean_Square_Error_Cal(equalization__chtaps,equalization)
MSE_T(k)=Mean_Square_Error_Cal(equalization__chtaps,equalization2)
%% Multiple Antennas(SIMO) part2
equalization_total=[equalization2 equlization_RX2 equlization_RX3 equlization_RX4];
equalization_mean=mean(equalization_total,2);
signal_rx_est_freq=((signal_rx_fx)./(equalization_mean));
signal_rx_fx=signal_rx_est_freq;


% delete the preamble
signal_rx_fx=[signal_rx_fx(:,3:end)];
%% CFO tracking with Pilot 
for i=1:data_L
mean_phase=Compute_Phase_CFO(signal_rx_fx(:,i),N_subcrr);%30*1
signal_rx__fine_CFO(:,i)=signal_rx_fx(:,i)*conj(mean_phase);
end
% Parallel to serial
[p, q] = size(signal_rx_fx);
mul_rx = p*q;
symb_rx = reshape (signal_rx_fx,mul_rx, 1);
% Demodulation
     [bits_rx] = demapping(symb_rx,Nbps);
 disp('Demapping done')

%% Results
% Bit Error Rate:
bitErrorRate = sum(abs(bits_tx - bits_rx),'all') / length(bits_tx);
disp('$$ Displaying results:');
disp(['BER:', num2str(bitErrorRate)]);
BER_total(k)=bitErrorRate;

%%  We display the constalation Tx and Rx
figure;
subplot(1,2,1); plot(real(symb_tx),imag(symb_tx),'rx'); 
title('Tx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
subplot(1,2,2); plot(real(symb_rx),imag(symb_rx),'g.'); 
title('Rx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])

end
close all;
figure;
semilogy(BER_total,'-x');
ylim([10^(-4) 10^(0)]);
hold on;
title('Bit error rate by the channel estimation in time-dom');
% title(['Bit error rate with length of cp is ',num2str(Lcp),' & data is ',num2str(data_L)]);
grid on;
xlabel('Varying SNR');
ylabel('Bit error rate');

%% MSE figure for channel estimation accuracy as function of SNR
figure;
plot(MSE_F,'-x');
hold on 
plot(MSE_T,'-x');
title('MSE comparision');ylabel('MSE');xlabel('list for a varying SNR');
legend('MSE in freq domain','MSE in time domain')



%% function 
function [I,CFO_est]=auto_corr(signal_rx)
% return the estimated Time acquisition and Frequency acquisition
[~, length2] =size(signal_rx);
num=1088;
T=1/(100e6);
% compute the auto-corrlation
for i=1:(length2-num-num)
temp(i)=sum(conj(signal_rx(i:i+num-1)).*(signal_rx(i+num:i+num+num-1)),2);
end

[M,I]=max(temp);
CFO_est=angle(temp(I))*(-1)/num/T;
%I is the start number of premble
figure;
plot(abs(temp))
xlabel('starting number of OFDM symbol ');ylabel('auto-corrlation value')
end

%%
function result=add_STO(signal_tx)
num_zero=12;
zero_delay=zeros(1,num_zero);
result=[zero_delay,signal_tx(1:end)];
%plot(result)
end

%%
function result=add_CFO(signal)
cfo_ppm = 12.75;
f_dc       = 3.6e9;
B          = 100e6;    
cfo_fr = (cfo_ppm/1e6)*f_dc;% ppm (part per million)
delta_omega=2*pi*cfo_fr;
T=1/(B);
for i=1:length(signal)
    signal_rx_update_receive_CFO(i)=signal(i)*exp(1j*delta_omega*i*T);
end
result=signal_rx_update_receive_CFO;
end

%%
function result=Compute_Phase_CFO(signal,num_sub)
sum=0;
for column=1:(num_sub/8):num_sub
sum=signal(column)+sum;
end
result=sum/8;
end

%%
function result=Multipath_Propa(signal_tx,chtaps)
signal_tx=conv(signal_tx,chtaps);
result=signal_tx(1:end-length(chtaps)+1);
end
%% 
function result=End_Delete_CP(signal_rx,m,n,length)
S2P_rx = reshape(signal_rx, m, n);
signal_rx_time_RCP = S2P_rx (length + 1:end,:);
result = fft(signal_rx_time_RCP,[],1);
end

%%
function result=Mean_Square_Error_Cal(tx,rx)
% calculate the MSE by the tx and rx
[column,~]=size(tx);
mse=0;
for i=1:column
    mse=abs((tx(i,1)-rx(i,1))^2)+mse;
end
result=1/column*mse;
end

%%
function result=compensating(signal,CFO_est,T,time_acquisition,length_constant)
signal_rx=signal;

% adjust the length of signal_rx to ensure no error during reshape later
if  (length(signal_rx)-time_acquisition+1)>=length_constant
    signal_rx=signal_rx(1+length(signal_rx)-length_constant:end);
else
    signal_rx=signal_rx((end-length_constant+1):end);
end

% compensating the correspond CFO estimation
    for j=1:length(signal_rx)
        result(j)=signal_rx(j)*exp(1i*CFO_est*j*T);
    end
end

%%
function [equalization2, time_acquisition, CFO_est]=SIMO(signal_tx,chtaps,preamble_tx)
% Multipath
signal_tx=Multipath_Propa(signal_tx,chtaps);
length_constant=length(signal_tx);


B=100e6;
Lcp=64;
N_subcrr=1024;
SNR_list = [20];
data_L=30;
Nbps=2;
T=1/(B);
Nbits = data_L * N_subcrr * Nbps;
EbN0 = SNR_list;% SNR
energy_of_signal=1/B * trapz(abs(signal_tx).^2);
Eb = energy_of_signal / Nbits ;
N0 = Eb/10.^(EbN0/10);
NoisePower = 2*N0*B;
signal_STO=add_STO(signal_tx);
noise = sqrt(NoisePower/2)*(randn(length(signal_STO),1)+...
1i*randn(length(signal_STO),1))';


%STO % CTO part
signal_STO=add_STO(signal_tx);
signal_STO_noise=signal_STO+noise;
signal_rx_STO_CFO=add_CFO(signal_STO_noise);
[time_acquisition, CFO_est]=auto_corr(signal_rx_STO_CFO);
signal_rx=compensating(signal_rx_STO_CFO,CFO_est,T,time_acquisition,length_constant);
m=N_subcrr+Lcp;n=data_L+2;
signal_rx_fx=End_Delete_CP(signal_rx,m,n,Lcp);
preamble_rx1 = signal_rx_fx(:,1);
preamble_rx2 = signal_rx_fx(:,2);
preamble_rx = [preamble_rx1 preamble_rx2];
CH_F = mean(preamble_rx./preamble_tx,2);
CH_T = ifft(CH_F);
CH_T = [CH_T(1:64);zeros((N_subcrr-Lcp),1)];
CH_F_t = fft(CH_T);


% CH_F(Channel estimation in Freq domain)
equalization = CH_F;
signal_rx_est_freq=((signal_rx_fx)./(equalization));
signal_rx_fx=signal_rx_est_freq;


% CH_T
equalization2= CH_F_t;
% signal_rx_est_time=((signal_rx_fx)./(equalization2));
% signal_rx_fx=signal_rx_est_time;
% delete the preamble
signal_rx_fx=[signal_rx_fx(:,3:end)];

% CFO tracking with Pilot 
for i=1:data_L
mean_phase=Compute_Phase_CFO(signal_rx_fx(:,i),N_subcrr);%30*1
signal_rx__fine_CFO(:,i)=signal_rx_fx(:,i)*conj(mean_phase);
end


%% Parallel to serial
% [p, q] = size(signal_rx__fine_CFO);
[p, q] = size(signal_rx_fx);
mul_rx = p*q;
symb_rx = reshape (signal_rx_fx,mul_rx, 1);

[bits_rx] = demapping(symb_rx,Nbps);
result=bits_rx;
disp('Demapping done')


end


