function distortedarray = add_CFO(inputarray, CFO, Ts)
%ADD_CFO Summary of this function goes here
%   Detailed explanation goes here
    distortedarray = zeros(size(inputarray));
    t = Ts*[0:length(inputarray)-1];
    distortedarray = inputarray.*exp(1j*2*pi*CFO*t);
end

