%% Modified methaneModel1Box_mrn, based on methaneModel1Box_mrn1.m by Mindy
% take out bio, ag, ff. 1box doesn't include boreal/nonboreal ratio of
% burning. Amount of geologic as pre-industrial methane source is debated (Thornton et al. 2021) 

ch4=importdata('ch4_adddata.mat'); % epica edml dome c, global methane
d13=importdata('EDML_d13CH4_20-80kyrBP_cut.txt'); % epica edml dome c, global isotopic ratio of methane
c2h2=importdata('data_4REU_20230710.mat'); % acetylene data measured by Jenn and Murat from SPICEcore and GISP2
o18=importdata('gkinis2020_cut.prn'); % from RECAP greenland ice core project
sorted_c2h2_ant1=importdata('sorted_c2h2_ant1.txt'); % additional acetylene measurements from jenn
md=importdata('c2h2_g2_md.txt'); % mid-depth of additional DO8 acetylene data
all_act_grn=importdata('all_act_grn.txt'); % old + new DO8 acetylene data from Jenn
% c2h2_add=importdata('Melter GISP2D 14Oct2022_LOG_QSC_20230801_ice.xls');
% d_ant_1=importdata('d_ant_1.txt'); %from jenn
% d_grn_1=importdata('d_grn_1.txt'); %from jenn
% o18_gkinis2020_b2k_bp_o18=importdata('o18_gkinis2020_b2k_bp_o18.txt');
% act_ppt_g2=importdata('act_ppt_g2.txt'); % new acetylene DO8 data

% sort acetylene data by age
% sorted_c2h2_grn=sortrows(d_grn_1); 
% sorted_c2h2_ant=sortrows(d_ant_1);
% sorted_c2h2_grn2=sortrows(act_ppt_g2);

% source histories in Tg/yr CH4 calling S from methaneModel1Box_mrn.m
S.year =  [35760 35895 36190 36485 36780 37075 37370 37665 37960 38255 38550]; % 295 intervals
S_pts=35500:295:38500;

grn_bb = 106.8; % (grn acetylene sensitivity from tracer model (ppt/Tgy^-1) (Nicewonger et al 2020), ppt/Tg/yr)
boreal_act=0.3; % grn - boreal dry matter burned emission factor g/kg to convert from acetylene (Nicewonger et al 2020)
%boreal_ch4=5.5; % grn

ant_bb = 6.8; % (antarctica acetylene sensitivity (ppt/Tgy^-1) (Nicewonger et al 2020), ppt/Tg/yr)
global_act = 0.33; % (global dry matter burned emission factor g/kg (Nicewonger et al 2020))
global_ch4 = 3.70; % (g/kg acetylene sensitivity for antarctica equation only)

mol_atmo = 1.77*10^20; % moles of "air" in the atmosphere (air = 29g/mol)
L_OH = 9.5; % methane lifetime (yrs) 

%% input of CH4 emissions
% Q: when simulating just BB, why is isotopic signature showing up at -17 and not -22.3?
% A: steady state equation: dx/dt = O = P - L -> E - k*x
% -> alpha = ^13K/^12K (fractionation rate)
% when inputting only bb, loss for 13C is a bit heavier than 12C, so it's 
% lost slower than 12C. This is why the model outputs bb a little heavier 
% than what is defined plainly as the isotopic end member for bb. 
% Because it includes k OH and k alpha (the fractionation rate), the 13C
% ratio ends up a little heavier (around -17) than the defined -22.3 in
% Mindy's 1box model

% using only bb and micro sources to optimize methane emissions
S_calc.bb = [ 43.7 44 46.5 49.5 53 58.5 55.3 53 47 42.5 43 ];
S_calc.micro = [ 102.3 105 107.5 114.5 127 133.5 134.8 131 118 107.5 108 ];
S_calc.geo = [ 0 0 0 0 0 0 0 0 0 0 0 ];

% keeping geo sources consistent within glacial values from Bock et al 2017
% S_calc.bb = [ 10 0.1 13 15 10 8 8 8 8 8 9 ]; %
% S_calc.micro = [ 83 98.9 89 97 123 127 137 128 106 90 92 ];
% S_calc.geo = [ 50 50 50 50 50 50 50 50 50 50 50];

% used acetylene DM data to optimize the methane bb and DM, then kept 
% the methane bb constant and changed geo and
% micro to fit delt 13C. does that work? it fits
% within the Bock 2017 literature min and max values for these sources:
% S_calc.bb = [ 10 0.1 13 15 11 8 9 8 7 10 7 ]; %
% S_calc.micro = [ 64 56.9 72 79 83 82 84 83 78 72 74 ];
% S_calc.geo = [ 65 90 67 68 88 95 95 95 84 66 70];

S_calc.ff = [ 0 0 0 0 0 0 0 0 0 0 0 ];
S_calc.bio = [ 0 0 0 0 0 0 0 0 0 0 0 ];
S_calc.ag = [ 0 0 0 0 0 0 0 0 0 0 0 ];

%% to loop through only 'new' literature value for 13C source isotope ratios
% calling from methaneModel1Box_mrn.m
I.iso_flag = {};
for i=1:11 % i = time steps (years)
    I.iso_flag(i) = {'new'};
end

for i = 1:11 % i = time steps for each source (years)
    S.geo=S_calc.geo(i); %the value of S for i-ith step
    S.micro=S_calc.micro(i);
    S.bb=S_calc.bb(i);
    S.ff=S_calc.ff(i);
    S.ag=S_calc.ag(i);
    S.bio=S_calc.bio(i);
        [CH4_Mtot, CH4_delC, CH4_delD] = methaneModel1Box_mrn1(S, L_OH, I.iso_flag); %run CH4 box model
        % Outputs from model for i-th year
        Out.CH4_ppb = (CH4_Mtot*10^12/mol_atmo)*1e9;
        Out.CH4_delC = CH4_delC;
        Out.CH4_delD = CH4_delD;
        % Variable for all runs
            results.CH4_ppb(:,i) = Out.CH4_ppb; 
            results.CH4_delC(:,i) = Out.CH4_delC; 
end

%% Interpolate to get GA for new acetylene points
%interp1(mid depth, ga from ch4, depth_pts of acetylene data points)

depth_pts = md; %number of depth points to iterate through interp, using mid depth (TD+BD/2)

ga_act=interp1(c2h2.methane.g2.depth, c2h2.methane.g2.gasage, depth_pts, 'linear');
S_bb_interp=interp1(S.year, S_calc.bb, S_pts, 'linear'); %output fitted curve from input bb emissions
S_micro_interp=interp1(S.year, S_calc.micro, S_pts, 'linear');
S_geo_interp=interp1(S.year, S_calc.geo, S_pts, 'linear');
% c2h2_interp_grn=interp1(sorted_c2h2_grn(:,1)*1000, sorted_c2h2_grn(:,2),S_pts, 'linear'); old grn data
c2h2_interp_grn_full=interp1(all_act_grn(:,1)*1000, all_act_grn(:,2), S_pts, 'linear'); % old and new additional grn data
c2h2_interp_ant=interp1(sorted_c2h2_ant1(:,1)*1000,sorted_c2h2_ant1(:,2), S_pts, 'linear');

%% EQ to convert acetylene (ppt) from Mindy's paper to biomass burning, then convert bb to methane emitted (tg) from biomass burning
% conv_ch4bb_grn = (c2h2_interp_grn_full)*(1/grn_bb)*(1e12/1)*(1/boreal_act)*(boreal_ch4)*(1/1e12);
% conv_ch4bb_ant = (c2h2_interp_ant)*(1/ant_bb)*(1e12/1)*(1/global_act)*(global_ch4)*(1/1e12); 
% global_conv_ch4bb = conv_ch4bb_grn+conv_ch4bb_ant;

% dry matter burned
conv_act_dm_grn = (c2h2_interp_grn_full)*(1/grn_bb)*(1e12/1)*(1/boreal_act)*(1/1e9); 
conv_act_dm_ant = (c2h2_interp_ant)*(1/ant_bb)*(1e12/1)*(1/global_act)*(1/1e9);
conv_ch41box_dm = ((S_bb_interp)*(1e12/1)*(1/global_ch4)*(1/1e9));
global_conv_act_dm = conv_act_dm_grn+conv_act_dm_ant;


%% Colors to plot

htmlGray = [178 178 178]/255;
darkblue = [0 153 167]/255;

%% Plot 

% Measured DO8 Acetylene (Jenn)
figure(1);
% subplot(2,2,1);
% plot(all_act_grn(:,1), all_act_grn(:,2), '.-r', 'linewidth', 0.75, 'MarkerSize',12)
% 
% hold on
% plot(sorted_c2h2_ant(:,1), sorted_c2h2_ant(:,2),  '.-k','linewidth', 0.75, 'MarkerSize',12)
% hold off
% 
% xlim([35.7 38.5])
% ylim([-5 70])
% title('Acetylene Concentration during DO8')
% ylabel('Acetylene (ppt)')
% xlabel('Ka Year (B.P.)')
% grid on
% set(gca, 'XDir','reverse')
% legend('Greenland', 'Antarctica')

%% Plot Literature CH4 Concentration, Optimized CH4 Emissions (EPICA EDML) and O18 G Kinis 2020
subplot(2,2,1)

yyaxis('right')
plot(o18(:,1)*1e3,o18(:,2), 'Color', htmlGray, 'linewidth', 0.55)
xlim([3.57*1e4 3.85*1e4]) 
ylabel('δ^1^8O (‰)')

yyaxis('left')
plot(S.year, results.CH4_ppb, '.-b', 'MarkerSize',12); % optimized CH4
xlim([3.57*1e4 3.85*1e4])
ylim([400 650])
grid on

hold on
plot(ch4.ga,ch4.ch4, '.', 'Color', darkblue, 'MarkerSize', 12) % literature values
hold off

title('Global CH_4 Total During DO8')
ylabel('CH_4 (ppb)')
xlabel('Year B.P.')
set(gca, 'XDir','reverse')
legend('Optimized CH_4', 'CH_4', 'δ^1^8O')

%% Plot Optimized Source Emissions of CH_4 from BB and Wetlands
subplot(2,2,2)

ylabel('Emissions (Tg of CH_4)')
ylim([0 200])
xlim([35760 38550])
title('Optimized Source Emissions of CH_4 from BB and Wetlands')
plot(S.year, S_bb_interp, '.-r', 'MarkerSize',12); % interpolated BB emissions
grid on

hold on
ylim([0 200])
title('Optimized Source Emissions of CH_4 from BB and Wetlands')
xlim([35760 38550])
ylabel('Emissions (Tg of CH_4)')
xlabel('Year B.P.')
plot(S.year, S_micro_interp, '.-b', 'MarkerSize', 12) %interpolated micro emissions
plot(S.year, S_geo_interp, '.-g', 'MarkerSize', 12)
hold off
legend('Biomass Burning','Wetlands', 'Geologic', Location='east');
set(gca, 'XDir','reverse')

%% Plot Isotopic Methane and Optimized Isotopic Methane

% yyaxis('right')
subplot(2,2,3)
plot(d13(:,2)*1e3,d13(:,3), '.r','linewidth', 0.75, 'MarkerSize',12) % literature iso data
ylim([-47 -44])
grid on

hold on
plot(S.year, results.CH4_delC, '.-b', 'MarkerSize',12); % optimized isotope
hold off

xlim([3.57*1e4 3.85*1e4])
xlabel('Year B.P.')
ylabel('δ^1^3C CH_4')
set(gca, 'XDir','reverse')
title('δ^1^3C CH_4')
legend('δ^1^3C', 'Optimized δ^1^3C', Location='southeast')

%% Plot Dry Matter Burned from Acetylene and Methane fire proxies

subplot(2,2,4)
% jenn/mindy's acetylene sensitivities based ch4 emissions
plot(S_pts,conv_ch41box_dm*1e-3, '.-k','linewidth', 0.75, 'MarkerSize',12) % 1e-3 to convert Tg to Pg (petagram)
ylim([0 18])
xlim([3.57*1e4 3.85*1e4])
grid on
ylabel('DM Burned (Pg/year)')
xlabel('Year B.P.')
title('Total Dry Matter (DM) Burned')

hold on
plot(S_pts, global_conv_act_dm*1e-3, '.-b','linewidth', 0.75, 'MarkerSize',12);
plot(S_pts, conv_act_dm_ant*1e-3, '.-r','linewidth', 0.75, 'MarkerSize',12);
plot(S_pts, conv_act_dm_grn*1e-3, '.-g', 'linewidth', 0.75, 'MarkerSize',12);
hold off

% hold on
% yyaxis('right') % global 1box methane model BB emissions
% plot(S_pts, S_calc.bb, '.-k','linewidth', 0.75, 'MarkerSize',12)
% ylabel('Act. conv. CH_4 emissions (Tg/year)')
% hold off
set(gca, 'XDir','reverse')
legend('CH_4-Inferred', 'Global Act-Inferred', 'Nonboreal', 'Boreal', Location='east')

%% Plot Acetylene-derived Biomass Burning in Tg of CH4 against Acetylene ppt

% subplot(4,1,4)
% yyaxis('left')
% plot(sorted_c2h2_grn(:,1)*1e3,conv_ch4_grn, '.-c','linewidth', 0.75, 'MarkerSize',12)
% xlim([3.57*1e4 3.85*1e4])
% set(gca, 'XDir','reverse')
% ylim([20 70])
% grid on
% ylabel('CH4 of bb (tg)')
% xlabel('ka Year B.P.')
% title('C_2H_2 (ppt) & Converted CH4 (Tg) from bb')
% 
% hold on
% yyaxis('right')
% plot(sorted_c2h2_grn(:,1)*1e3,sorted_c2h2_grn(:,2), '.-b','linewidth', 0.75, 'MarkerSize',12)
% xlim([3.57*1e4 3.85*1e4])
% set(gca, 'XDir','reverse')
% grid on
% ylabel('C_2H_2 GRN (ppt)')
% ylim([20 70])
% xlabel('ka Year B.P.')
% title('C_2H_2 (ppt) & Converted CH4 (Tg) from bb')
% hold off
