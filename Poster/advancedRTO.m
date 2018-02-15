
% generate advanced RTO pictures
% starting with transient

% A+B->C
% r = k*cA*cB

syms F_A F_B F_C

model.k = 5;
plant.k = 6;

phi = 10*F_C - 3*F_B - 2*F_A;







