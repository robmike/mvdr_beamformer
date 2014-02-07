function ber = get_ber_approx(energy,coeffs,rakevec,corrmat)

% Approximation to bit error rate using SINR
% corrmat - Correlation matrix of interference + noise

ber = energy*abs(coeffs'*rakevec)^2 / (coeffs'*corrmat*coeffs);
ber = real(ber); % Computations introduce imaginary components
ber = Q(ber);