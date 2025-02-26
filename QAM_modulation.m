function symbols = QAM_modulation(bits,bitspersymbol)
%QAM_MODULATION Heavily inspired 
%   Detailed explanation goes here

    am_bits = size(bits,1);
    if mod(am_bits,bitspersymbol) ~= 0
        symbols = 0;
        disp('error: bit padding has not been implemented as of yet')
        return
    end
    
    am_symb = am_bits/bitspersymbol;
%     if mod(am_symb,2) ~= 0
%         symbols = 0;
%         disp('error: symbol padding has not been implemented as of yet')
%         return
%     end
    
    if mod(bitspersymbol,2) ~= 0
        symbols = 0;
        disp('error: uneven QAM design not implemented yet')
        return
    end
    re_bitspersymbol = bitspersymbol/2;
    im_bitspersymbol = re_bitspersymbol;
    %split bits into two
    re_bits = bits(1:am_bits/2);
    im_bits = bits(am_bits/2 +1: end);

  
    re_int = bit2int(re_bits, re_bitspersymbol);
    im_int = bit2int(im_bits, im_bitspersymbol);
    
    %and now to a QAM symbol (calculations taken from Prof. Horlin's
    %Communication Channel slides
    re_sig = sqrt( sum(([0:2^re_bitspersymbol-1]-(2^re_bitspersymbol-1)/2).^2) /2^re_bitspersymbol ); %symbol spreading factor
    im_sig = sqrt( sum(([0:2^im_bitspersymbol-1]-(2^im_bitspersymbol-1)/2).^2) /2^im_bitspersymbol ); %symbol spreading factor
    re_symbols = 1/(re_sig*sqrt(2)) * (re_int - (2^re_bitspersymbol-1)/2); % - (2^re_bitspersymbol-1)/2 to center it around origin
    im_symbols = 1/(im_sig*sqrt(2)) * (im_int - (2^im_bitspersymbol-1)/2); % - (2^re_bitspersymbol-1)/2 to center it around origin
    symbols = re_symbols + 1j*im_symbols;

end

