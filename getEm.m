
function [eint] = getEm(tseg, e, tinterp)
ee = [e(1) e(2) e(3) e(4) e(5) e(6) e(7) e(8) e(9) e(10) e(11) e(12)];
eint = interp1( tseg, ee, tinterp, 'linear', 'extrap');  
end

%***************************************************
%***************************************************
% Calculates lineraly interpolated emissions for gas species
% Inputs 
%   tseg time segments [1000 1500 1700 1800 1900 2000]
%   e is emissions for [1000-1500, 1700-1800, 1900, 2000] is a 4x1
% Outputs
% eint is linearly interpolated emisisons over tinterp
%***************************************************
%***************************************************
%Written by Mindy & Eric, 08Feb2018