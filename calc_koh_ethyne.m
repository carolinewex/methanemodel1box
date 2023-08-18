function [box_kOH] = calc_koh_ethyne(lat_top, lat_bot)
%************************************************************
%calculates k_OH for acetylene 
%inputs: 
%   lat_top is nX1 vector in degrees top of each box
%   lat_bot is nX1 vector in degrees bot of each box 
%outputs:
%   box_kOH is nX1 vector of loss frequencies in units yr^-1
%*************************************************************

load('Spiv2000.mat') 
p_Pa = p.*100; %converting mb to Pa    
R = 8.3144598*10^6; %gas constant cm^3*Pa/K*mol
avo = 6.022*10^23; %avogadro's #  

molec_cm3 = zeros(9,24); 

for i = 1:length(p_Pa)
    for j = 1:length(lat)
    
    molec_cm3(i,j) = (p_Pa(i)/(R*globalT(i,j)))*avo; 
    end
end

molec_cm3(isinf(molec_cm3)) = 0; 

k0 = (5.5*10^-30); %k0(T) = ko*(T/300)^-0 (so no temp dependence)
kinf = (8.3*10^-13).*((globalT./300).^2); %kinf(T) = kinf*(T/300)^2

a = (k0.*molec_cm3); 
b = (a./kinf);
%c = (k0.*molec_cm3)./(kinf.*(globalT./300).^2);
c = log10(b).^2;

k_eff = 1.05*(a./(1+b)).*0.6.^(1./(1+c));

% plot(k_eff); 

% k_low = a; 
% k_inf = 
% global_k_OH_acetylene = globalOH.*1e5.*k_low;
% global_k_OH_acetylene = globalOH.*1e5.

global_k_OH_acetylene = globalOH.*1e5.*k_eff;

%plot(global_k_OH_acetylene)

 %global_k_OH_acetylene = globalOH.*1e5*(4.15*10^-13).*((globalT./300).^-2); %high pressure limi
%global_k_OH_acetylene = globalOH.*1e5*(8.3*10^-13).*((globalT./300).^2); %high pressure limit
%   from JPL - 
%global_k_OH_acetylene = globalOH.*1e5.*(9.4*10^-12).*exp(-700./globalT); %Atkinson et al 1990 

  
 % ko(T) = A'*(T/300)^-n in cm3 molecule^-1 sec^-1
 %global_k_OH_ethane = globalOH.*1e5.*8.7e-12.*exp(-1070./globalT);  %from Sander et al. 2006, JPL publication pg. 1-11 cm3 molecule-1 sec-1
 
    %K value from JPL evaluation #17, 2011
    p_frac = [-diff(p); 0]./sum(-diff(p)); 
    p_frac2 = repmat(p_frac, 1, size(globalOH, 2)); 
    k_zonal = nansum(global_k_OH_acetylene.*p_frac2, 1)'; 
    SA = 4*pi*(6378100^2);
    lat(1) = -90;
    lat(end)= 90;
    lat_i = (-90:1:90)';
    k_zonal_i = interp1(lat, k_zonal, lat_i);
    area_i = cosd(lat_i);
    box_kOH = []; 
    for i = 1:length(lat_top)
        vr= find(lat_i >= lat_bot(i) & lat_i <=lat_top(i)); %indices of latitudes 90-60N
        box_kOH(i) = sum(k_zonal_i(vr).*area_i(vr))./sum(area_i(vr)); %rate constant s-1
    end
    box_kOH = box_kOH * 60*60*24*365;
end 