function ToA_est_bits = generate_ToA_est_bits(symbol_length)
    %idea: alternating ones and zeroes, to xcorr
    ToA_est_bits = zeros(symbol_length,1);
    index = 1;
    addition_term = 1;
    while index <= symbol_length
        ToA_est_bits(index) = 1;
        addition_term = addition_term + 1;
        index = index + addition_term;
    end
%     for i = 2:2:symbol_length
%         ToA_est_bits(i) = 1; %% more compute-efficient way?
%     end  
end