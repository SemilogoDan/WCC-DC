function [Symbols_fixed, total_indices] = add_pilots(Symbols, Q, pilot_symbol, pilot_percentage)
%ADD_PILOTS Summary of this function goes here
%   Detailed explanation goes here
Symbols_fixed = Symbols;
pilot_amount_per_row  = floor(Q*pilot_percentage);
am_rows = floor(length(Symbols)/Q);
pilot_indices_row = floor(linspace(1,ceil(Q*(pilot_amount_per_row-1)/pilot_amount_per_row), pilot_amount_per_row));
total_indices = zeros(pilot_amount_per_row*am_rows,1);
for i =0:am_rows-1
    Symbols_fixed(pilot_indices_row + i*Q) = pilot_symbol;
    total_indices(i*pilot_amount_per_row+1:(i+1)*pilot_amount_per_row) = pilot_indices_row + i*Q;
end
end