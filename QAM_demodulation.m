% function bits = QAM_demodulation(symbols,bitspersymbol)
% %QAM_DEMODULATION Summary of this function goes here
% %   Detailed explanation goes here
%     re_symbols = real(symbols);
%     im_symbols = imag(symbols);
%     re_bitspersymbol = bitspersymbol/2;
%     im_bitspersymbol = re_bitspersymbol;
% 
%     %and now to a integers (calculations taken from Prof. Horlin's
%     %Communication Channel slides
%     re_sig = sqrt( sum(([0:2^re_bitspersymbol-1]-(2^re_bitspersymbol-1)/2).^2) /2^re_bitspersymbol ); %symbol spreading factor
%     im_sig = sqrt( sum(([0:2^im_bitspersymbol-1]-(2^im_bitspersymbol-1)/2).^2) /2^im_bitspersymbol ); %symbol spreading factor
%     re_int = re_sig * sqrt(2) * re_symbols + (2^re_bitspersymbol-1)/2;
%     im_int = im_sig * sqrt(2) * im_symbols + (2^im_bitspersymbol-1)/2;
% 
%     % the received values of course are not intact due to interference
%     re_int_est = round(re_int);
%     im_int_est = round(im_int);
% 
%     %remove bogus results
%     for i = 1:length(re_int_est)
%         if re_int_est(i) < 0
%             re_int_est(i) = 0;
%         elseif re_int_est(i) > 2^re_bitspersymbol - 1
%             re_int_est(i) = 2^re_bitspersymbol - 1;
%         end
%     end
%     for i = 1:length(im_int_est)
%         if im_int_est(i) < 0
%             im_int_est(i) = 0;
%         elseif im_int_est(i) > 2^im_bitspersymbol - 1
%             im_int_est(i) = 2^im_bitspersymbol - 1;
%         end
%     end
%     re_bit = int2bit(re_int_est, re_bitspersymbol);
%     im_bit = int2bit(im_int_est, im_bitspersymbol);
% 
%     bits = [re_bit,im_bit];
% end
% 
% 
% 

function bits = QAM_demodulation(symbols, bitspersymbol)
    re_symbols = real(symbols);
    im_symbols = imag(symbols);
    re_bitspersymbol = bitspersymbol/2;
    im_bitspersymbol = re_bitspersymbol;

    % and now to integers (calculations taken from Prof. Horlin's Communication Channel slides
    re_sig = sqrt( sum(([0:2^re_bitspersymbol-1]-(2^re_bitspersymbol-1)/2).^2) /2^re_bitspersymbol ); % symbol spreading factor
    im_sig = sqrt( sum(([0:2^im_bitspersymbol-1]-(2^im_bitspersymbol-1)/2).^2) /2^im_bitspersymbol ); % symbol spreading factor
    re_int = re_sig * sqrt(2) * re_symbols + (2^re_bitspersymbol-1)/2;
    im_int = im_sig * sqrt(2) * im_symbols + (2^im_bitspersymbol-1)/2;

    % the received values of course are not intact due to interference
    re_int_est = round(re_int);
    im_int_est = round(im_int);

    % remove bogus results
    re_int_est(re_int_est < 0) = 0;
    re_int_est(re_int_est > 2^re_bitspersymbol - 1) = 2^re_bitspersymbol - 1;
    im_int_est(im_int_est < 0) = 0;
    im_int_est(im_int_est > 2^im_bitspersymbol - 1) = 2^im_bitspersymbol - 1;

    % Truncate or pad arrays for size compatibility
    if length(re_int_est) > length(re_symbols)
        re_int_est = re_int_est(1:length(re_symbols));
    elseif length(re_int_est) < length(re_symbols)
        re_int_est = [re_int_est; zeros(length(re_symbols) - length(re_int_est), 1)];
    end
    
    if length(im_int_est) > length(im_symbols)
        im_int_est = im_int_est(1:length(im_symbols));
    elseif length(im_int_est) < length(im_symbols)
        im_int_est = [im_int_est; zeros(length(im_symbols) - length(im_int_est), 1)];
    end

    re_bit = int2bit(re_int_est, re_bitspersymbol);
    im_bit = int2bit(im_int_est, im_bitspersymbol);

    bits = [re_bit; im_bit];
end
