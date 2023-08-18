function [Out] = calc_emissions(I)

load('S.mat'); %this is the RCP emissions inventory , how do we 'turn off' rcp scenarios when its connected to a bunch of stuff in the code

dlat_box = 30; %width of each latitude box in model 
lat_top = [90:-dlat_box:-60]';
lat_bot = [60:-dlat_box:-90]';


[k_OH_ethane] = I.k_scaler*calc_koh_ethane(lat_top, lat_bot); 
[k_OH_ethyne] = I.k_scaler*calc_koh_ethyne(lat_top, lat_bot);%acetylene 
[k_OH_methane] = I.k_scaler*calc_koh_methane(lat_top, lat_bot); 

SA = 4*pi*(6378100^2); %total surface area of earth 
mol_atmo = 1.77*10^20; %moles of "air" in the atmosphere (air = 29g/mol)
mol_trop = 0.8*mol_atmo; %moles of "air" in the troposphere (up to 200 hPa, 80%)
mass_trop = mol_atmo*28.97; % mass of troposphere in grams (80% of mass) 

A1 = 2*pi*(6378100^2)*(sind(90)-sind(60)); %surface area in latitude bands 
A2 = 2*pi*(6378100^2)*(sind(60)-sind(30));
A3 = 2*pi*(6378100^2)*(sind(30)-sind(0));
areas = [A1/SA, A2/SA, A3/SA, A3/SA, A2/SA, A1/SA]; %ratio of SA in each lat band
scaler = [0.95, 0.95, 1.05, 1.05, 0.95, 0.95]; %mass in each box scaled based on trop. height (shallow in poles)
mass_box = mass_trop*areas.*scaler; % mass of air in each "box" [g] %To take into account the lower tropopause, polar and mid lat bands mass decreased by 5% while tropical bands increase 5%
moles_box = mol_trop*areas.*scaler;


%HERE IS THE FUNCTION THAT CALCULATES STEADY STATE MIXING
[A_ethane] = inversion_6box(k_OH_ethane); 
[A_ethyne] = inversion_6box(k_OH_ethyne);
[A_methane] = inversion_6box(k_OH_methane); 

C2H6_g_mol = 30.07; % C2H6 grams/mol 
C2H2_g_mol = 26.04;  
CH4_g_mol = 16.04; 

%Calculates global lifetime of methane based on the 6-box model & OH
%concentrations from Spivakosky, does the global lifetime of methane change
%based on what years it is? so are the numbers from Spivakosy irrelevant
%for 37000 years ago?
global_lftm_methane = 1/((mass_box(1)*k_OH_methane(1) + mass_box(2)*k_OH_methane(2) + mass_box(3)*k_OH_methane(3) + mass_box(4)*k_OH_methane(4) + mass_box(5)*k_OH_methane(5) + mass_box(6)*k_OH_methane(6))/sum(mass_box));

%Scales RCP fossil inventory to Tg/yr requested - turning of rcp doesn't
%matter as long as ff are turned off--bc that's the only thing being
%affected by rcp?
E.ff.act = S.RCP.act.emiss_tgyr_tot.*I.Tg_ff_act/sum(S.RCP.act.emiss_tgyr_tot)*10^12; 
E.ff.eth = S.RCP.eth.emiss_tgyr_tot.*I.Tg_ff_eth/sum(S.RCP.eth.emiss_tgyr_tot)*10^12; 
E.ff.CH4 = S.RCP.CH4.emiss_tgyr_tot.*I.Tg_ff_CH4/sum(S.RCP.CH4.emiss_tgyr_tot); 

%Biofuel inventory based on Yevick & Logan 2003 distribution of NMHCs
bio = [0.01, 0.309, 0.532, 0.142, 0.007, 0];

E.bio.act = I.Tg_bio_act.*bio*10^12;
E.bio.eth = I.Tg_bio_eth.*bio*10^12;
E.bio.CH4 = I.Tg_bio_CH4.*bio;

%Geologic emissions distribution based on Etiope and Ciccioli, 2009
geo = [0.05, 0.44, 0.33, 0.15, 0.03, 0];

E.geo.act = I.Tg_geo_act.*geo*10^12;
E.geo.eth = I.Tg_geo_eth.*geo*10^12;
E.geo.CH4 = I.Tg_geo_CH4.*geo;

%%%%GFED4%%%%

E.DM.ag = [63232952, 105637960000,	94001684000, 57478484000, 4350497300, 0];
E.DM.bor = [144705770000, 234652680000,	199660,	2571382, 0, 0];
E.DM.def = [0, 0, 185165380000,	401315590000,	0, 0];
E.DM.peat = [6462321200, 8308697600, 15514121000, 68251984000, 	0,	0];
E.DM.sav = [83993112, 54993875000, 1021639800000, 1823626100000, 16384378000, 0];
E.DM.tem = [244809536, 37365940000, 27763874000, 11989496000, 24871608000, 0];
E.DM.all = [E.DM.ag; E.DM.bor; E.DM.def; E.DM.peat; E.DM.sav; E.DM.tem]; 
E.DM.all_scaled = E.DM.all.*I.nonboreal_scale;
   
%%% biomass burning emission factors for each biome 
bb_ratio_eth = [0.91, 1.79, 0.71, 0.71, 0.66, 0.63]';
bb_ratio_act = [0.27, 0.18, 0.44, 0.06, 0.24, 0.26]';
bb_ratio_CH4 = 1.3*[5.82, 5.96, 5.07, 20.8, 1.94, 3.36]';

%making the bb. emission ratios the correct size matrix
bb_ratio_eth = repmat(bb_ratio_eth, 1, 6);
bb_ratio_act = repmat(bb_ratio_act, 1, 6);
bb_ratio_CH4 = repmat(bb_ratio_CH4, 1, 6);

emiss_bb_eth = E.DM.all.*bb_ratio_eth;
emiss_bb_act = E.DM.all.*bb_ratio_act;
emiss_bb_CH4 = E.DM.all.*bb_ratio_CH4;

%Totaling GFED emissions from diff. fire types for each "box"
for i = 1:6
   E.bb.eth(i) = sum(emiss_bb_eth(:,i)); 
   E.bb.act(i) = sum(emiss_bb_act(:,i)); 
   E.bb.CH4(i) = sum(emiss_bb_CH4(:,i))/10^12;     
end

boreal = [0.38 0.62 0 0 0 0]; 
%boreal fire contribution to each box
boreal_eth = boreal.*I.boreal_percent_eth*sum(E.bb.eth);
boreal_act = boreal.*I.boreal_percent_act*sum(E.bb.act); %based on GFED4s
boreal_CH4 = boreal.*I.boreal_percent_CH4*sum(E.bb.CH4); 

E.bb.eth_boreal = I.boreal_scale.*boreal_eth; 
E.bb.act_boreal = I.boreal_scale.*boreal_act; %scales boreal fires seperate from rest of fires
E.bb.CH4_boreal = I.boreal_scale.*boreal_CH4; 

E.bb.eth_nonboreal = I.nonboreal_scale.*(E.bb.eth-boreal_eth); 
E.bb.act_nonboreal = I.nonboreal_scale.*(E.bb.act-boreal_act);
E.bb.CH4_nonboreal = I.nonboreal_scale.*(E.bb.CH4-boreal_CH4);

%S.ag + S.geo + S.micro + S.bb; %THIS SENDS INFO TO THE METHANE MODEL! Need
%to sum up the 6 boxes of methane and send to the solver

%Inversion for ethane & acetylene
M_eth = A_ethane\(E.ff.eth + E.bio.eth + E.geo.eth + E.bb.eth + E.bb.eth_boreal)';
M_act = A_ethyne\(E.ff.act + E.bio.act + E.geo.act + E.bb.act + E.bb.act_boreal)';

%Inversion for methane!
E_CH4.ff = sum(E.ff.CH4);
E_CH4.geo = sum(E.geo.CH4); 
E_CH4.micro = sum(I.Tg_micro_CH4); 
E_CH4.ag = sum(I.Tg_ag_CH4); 
E_CH4.bb = sum(E.bb.CH4_nonboreal) + sum(E.bb.CH4_boreal); 
E_CH4.bio = sum(E.bio.CH4); 

L_OH = global_lftm_methane; 
[CH4_Mtot, CH4_delC, CH4_delD] = methaneModel1Box_mrn(E_CH4, L_OH, I.iso_flag);

%calulation of box mixing ratios for modern (M) and pre-industiral (M2)
Out.CH4_ppb = (CH4_Mtot*10^12/mol_atmo)*1e9;
Out.CH4_delC = CH4_delC;
Out.CH4_delD = CH4_delD;
Out.CH4_lftm = global_lftm_methane;
Out.budgets = E;
Out.eth_ppt = (1e12.*(M_eth./C2H6_g_mol)./moles_box');
Out.act_ppt = (1e12.*(M_act./C2H2_g_mol)./moles_box');
Out.M_eth = M_eth; 
Out.M_act = M_act; 
% Out.budgets_CH4 = E_CH4;
% Out.ethane_bb_nonboreal = E.bb.eth_nonboreal;
% Out.act_bb_nonboreal = E.bb.act_nonboreal;
% Out.ethane_bb_boreal = E.bb.eth_boreal; 
% Out.act_bb_boreal = E.bb.act_boreal; 
% Out.CH4_bb_boreal = E.bb.CH4_boreal; 
% Out.CH4_bb_nonboreal = E.bb.CH4_nonboreal; 

