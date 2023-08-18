function [box_kOH] = calc_koh(lat_top, lat_bot)
%************************************************************
%calculates k_OH for methane
%inputs: 
%   lat_top is nX1 vector in degrees top of each box
%   lat_bot is nX1 vector in degrees bot of each box 
%outputs:
%   box_kOH is nX1 vector of loss frequencies in units yr^-1
%
% 24 latitudes x 9 vertical layers 
%*************************************************************
    load('Spiv2000.mat') 
    

    %second order recation in Arrhenius form from JPL
    %k(T) = A*exp((-E/R)/T)
    
   global_k_OH_methane = globalOH.*1e5.*2.45e-12.*exp(-1775./globalT);  %from JPL
   %1e5 is a scaler for global OH gridded concentrations  
   
    p_frac = [-diff(p); 0]./sum(-diff(p)); %calculating the weighted pressure levels i.e. (1000mb - 900mb)/total mb of 900 = .11 
    p_frac2 = repmat(p_frac, 1, size(globalOH, 2)); %creates matrix with p_frac values in every latitude box (9x24) 
    k_zonal = nansum(global_k_OH_methane.*p_frac2, 1)'; %each vertical layer weighted and summed to get 24 zonal k_OH values 
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


%%modeled kOH in each box:
%90-60N = 1.4578 --> lifetime of 250 days
%60-30N = 4.1802 --> lifetime of 87 days
%30N-0 = 8.7496 --> lifetime of 41 days
%0-30S = 8.875 --> 41 days
%30-60S = 3.4869 --> 104 days
%60-90s = 1.038 --> 351 days 

%global avg. lifetime = 95 days (compared to Xiao et al. 2008 of 80 days)