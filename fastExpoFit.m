function T2s = fastExpoFit(m, deltaTE)
% implementation of Pei, "Algorithm for Fast Monoexponential Fitting Based on
% Auto-Regression on Linear Operations (ARLO) of Data", 2015
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25137
%
% Considers the corresponding erratum
% https://onlinelibrary.wiley.com/doi/10.1002/mrm.27807
%
% notation as in the publication
% m is the decaying signal
% deltaTE the spacing between echos


s = deltaTE/3*(m(1:(end-2)) + 4*m(2:(end-1)) + m(3:end));
delta = m(1:(end-2)) - m(3:end);

T2s = (sum(s.^2) + deltaTE/3*sum(s.*delta))/(deltaTE/3*sum(delta.^2) + sum(s.*delta));




