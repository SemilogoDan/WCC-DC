function preamble=form_premble()
Nbps=4;
% BPSK: generate a series of bits and use it for by premble by its BPSK  formation
bits = randi([0 1], [1 (N_subcrr)*4]);
bits_tx = bits' ;
preamble = mapping(bits_tx,Nbps,'qam');
save('preamble.mat','preamble');
end
