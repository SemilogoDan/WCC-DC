function ToA_est_signal_MIMO = generate_ToA_est_signal_received(symbol_length, cp_length, bitspersymbol, Ntx)
%GENERATE_TOA_EST_SIGNAL_RECEIVED Summary of this function goes here
%   Detailed explanation goes here
ToA_est_signal = generate_ToA_est_signal(symbol_length,cp_length, bitspersymbol);
ToA_est_signal_MIMO = repmat(ToA_est_signal, Ntx, 1);
end

