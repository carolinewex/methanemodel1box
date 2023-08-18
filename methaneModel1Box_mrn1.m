function [CH4_Mtot, CH4_delC, CH4_delD] = methaneModel1Box_mrn1(S, L_OH, iso_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: 
%S = methane source structure
%L_OH = methane lifetime (yrs)
%iso_flag = 'old' or 'new' for end-members
%OUTPUTS:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%methane model


global S13C S12C SD L

% %source histories in Tg/yr CH4
% S.year =  [ 1700 1800 1900 1920 1940 1960 1980 2000 ];
% S.ag =    [ 10 10 15 30 60 60 80 50 ];
% S.geo =   [ 10 10 19 30 60 115 125 100 ];
% S.micro = [ 220 220 220 220 240 275 275 275 ];
% S.bb =    [ 18 18 20 25 25 25 25 25 ];
S.totalTg = S.ag + S.geo + S.micro + S.bb + S.ff + S.bio;
S.totalTMol = S.totalTg/16; %approximate

%13C source isotope ratios
if strcmp('new', iso_flag)
    Rstd13C = 0.0112372; %Craig, 1957
    SdelC.ag = -62.3/1000; %uncertainity is %-60/1000;
    SdelC.geo = -43/1000;          %-40/1000;
    SdelC.micro =  -62.3/1000; %-61.7/1000;      %-60/1000;
    SdelC.bb = -22.3/1000; %-22.5/1000; %Shwietzke et al 2017       %-20/1000;
    SdelC.bio = -22.3/1000; %-22.5/1000; %Shwietzke et al 2017       %-20/1000;
    SdelC.ff = -43/1000; %Shwietzke et al 2017
else
    Rstd13C = 0.0112372; %Craig, 1957
    SdelC.ag = -60/1000; 
    SdelC.geo = -40/1000;         
    SdelC.micro =  -60/1000; 
    SdelC.bb = -20/1000; 
    SdelC.bio = -20/1000; 
    SdelC.ff = -40/1000;
end

S13Cratio.ag = (SdelC.ag + 1) * Rstd13C;
S13Cratio.geo =  (SdelC.geo + 1) * Rstd13C;
S13Cratio.micro = (SdelC.micro + 1) * Rstd13C;
S13Cratio.bb = (SdelC.bb + 1) * Rstd13C;
S13Cratio.ff = (SdelC.ff +1) * Rstd13C; 
S13Cratio.bio = (SdelC.bio +1) * Rstd13C; 

%check on source del 13C
sourceDel13C_1 = 1e3*(SdelC.ag*S.ag + SdelC.geo*S.geo + SdelC.micro*S.micro + SdelC.bb*S.bb + SdelC.bio*S.bio + SdelC.ff*S.ff)./16./S.totalTMol(1);
%compute isotope sources by solving simultaneous equations:
%   17*TMol13C + 16*TMol12 = TgTotal
%   TMol12C*R = TMol13C
S13C.agTMol = S.ag./(17+(16/S13Cratio.ag));
S13C.geoTMol = S.geo./(17+(16/S13Cratio.geo));
S13C.microTMol = S.micro./(17+(16/S13Cratio.micro));
S13C.bbTMol = S.bb./(17+(16/S13Cratio.bb));
S13C.ffTMol = S.ff./(17+(16/S13Cratio.ff));
S13C.bioTMol = S.bio./(17+(16/S13Cratio.bio));

S12C.agTMol = S.ag./((17*S13Cratio.ag+16));
S12C.geoTMol = S.geo./((17*S13Cratio.geo+16));
S12C.microTMol = S.micro./((17*S13Cratio.micro+16));
S12C.bbTMol = S.bb./((17*S13Cratio.bb+16));
S12C.ffTMol = S.ff./((17*S13Cratio.ff+16));
S12C.bioTMol = S.bio./((17*S13Cratio.bio+16));

S13C.totalTMol = S13C.agTMol + S13C.geoTMol + S13C.microTMol + S13C.bbTMol + S13C.ffTMol + S13C.bioTMol;
S12C.totalTMol = S12C.agTMol + S12C.geoTMol + S12C.microTMol + S12C.bbTMol + S12C.ffTMol + S12C.bioTMol;
%recheck source del 13C
%S12C.totalTMol = S.totalTMol - S13C.totalTMol;
sourceDel13C_2 = (((S13C.totalTMol./S12C.totalTMol)/Rstd13C)-1)*1000;

%13C source histories in Tg/yr 13CH4
S13C.ag = S13C.agTMol .* 17;
S13C.geo = S13C.geoTMol * 17;
S13C.micro = S13C.microTMol * 17;
S13C.bb = S13C.bbTMol * 17;
S13C.ff = S13C.ffTMol * 17;
S13C.bio = S13C.bioTMol * 17;
S13C.totalTg = S13C.ag + S13C.geo + S13C.micro + S13C.bb + S13C.ff + S13C.bio;

%12C source histories in Tg/yr 12CH4
S12C.ag = S12C.agTMol .* 16;
S12C.geo = S12C.geoTMol .* 16;
S12C.micro = S12C.microTMol  .* 16;
S12C.bb = S12C.bbTMol  .* 16;
S12C.ff = S12C.ffTMol  .* 16;
S12C.bio = S12C.bioTMol  .* 16;
S12C.totalTg = S12C.ag + S12C.geo + S12C.micro + S12C.bb + S12C.ff + S12C.bio;


%D source isotope ratios - ******CHECK THIS BEFORE USING
RstdD = 155.76 * 1E-6;  %D/H ratio in V-SMOW Gonfantini, 1978 should be increased by 0.2 per mil
SdelD.ag = -330/1000;
SdelD.geo = -175/1000;
SdelD.micro = -322/1000;
SdelD.bb = -169/1000;
SdelD.ff = -175/1000;
SdelD.bio = -169/1000;

SDratio.ag = (SdelD.ag + 1) * RstdD;
SDratio.geo = (SdelD.geo + 1) * RstdD;
SDratio.micro = (SdelD.micro + 1) * RstdD;
SDratio.bb = (SdelD.bb + 1) * RstdD;
SDratio.ff = (SdelD.ff + 1) * RstdD;
SDratio.bio = (SdelD.bio + 1) * RstdD;

%D source histories in Tg/yr D CH4
SD.ag = S12C.ag * SDratio.ag * (18/16);
SD.geo = S12C.geo * SDratio.geo * (18/16);
SD.micro = S12C.micro * SDratio.micro * (18/16);
SD.bb = S12C.bb * SDratio.bb * (18/16);
SD.ff = S12C.ff * SDratio.ff * (18/16);
SD.bio = S12C.bio * SDratio.bio * (18/16);

SD.totalTg = SD.ag + SD.geo + SD.micro + SD.bb + SD.ff + SD.bio;
SD.totalTMol = SD.totalTg/18;

%sinks - troposphere, stratosphere, soils
% L.year = [1700 1900 2000];
% %calculate OH scaling factor (assuming OH changed 25% due to 700-1700 ppm
% %increase
% TmolAtm = 5E18*1E3/(30*1E12); 
% L.scaling.slope = -0.25/((1.7-0.7)*1E-6*TmolAtm);
% L.scaling.intercept = 1 - (L.scaling.slope*0.7*1E-6*TmolAtm);

%L.scaling = [ 1 1 0.75 ]; now computed dynamically 
% L.C13alpha = 0.88*.9961 + 0.07*.9847 + .05*.9824; %13k/12k from Mischler et al GBC 2009
% L.k = 1/7.6; %years
% L.k13C = L.k * L.C13alpha; %note alpha refers to molar ratios
% L.Dalpha =  0.88*.7729 + 0.07*.7005 + .05*.8764; %Dk/Hk from Mischler et al GBC 2009
% L.kD = L.k * L.Dalpha;

% L.k.total = 1/7.6; %present day methane total lifetime
% L.k.OH = 0.88 * L.k.total;
% L.k.soil = 0.07 * L.k.total;
% L.k.strat = 0.05 * L.k.total;
% L.k.Cl = 0 * L.k.total; %zero Cl oxidation 

L.k.total = (1/L_OH)/.88; %present day methane total lifetime - Eric had (1/7.6)
L.k.OH = 0.88 * L.k.total;
L.k.soil = 0.07 * L.k.total;
L.k.strat = 0.05 * L.k.total;
L.k.Cl = 0 * L.k.total; %zero Cl oxidation 



L.C13alpha.OH = 0.9961;
L.C13alpha.soil = 0.9847;
L.C13alpha.strat = 0.9824; 
L.C13alpha.Cl = 1/1.0621; %Tyler et al. GRL 2000 at 298
L.k13C.OH = L.k.OH * L.C13alpha.OH;
L.k13C.soil = L.k.soil * L.C13alpha.soil;
L.k13C.strat = L.k.strat * L.C13alpha.strat;
L.k13C.Cl = L.k.Cl * L.C13alpha.Cl;

L.Dalpha.OH = 0.7729;
L.Dalpha.soil = 0.7005;
L.Dalpha.strat = 0.8764;
L.Dalpha.Cl = 1/1.474; %Tyler et al. GRL 2000 at 298
L.kD.OH = L.k.OH * L.Dalpha.OH;
L.kD.soil = L.k.soil * L.Dalpha.soil;
L.kD.strat = L.k.strat * L.Dalpha.strat;
L.kD.Cl = L.k.Cl * L.Dalpha.Cl;

%initial conditions in atmosphere - assume steady state
massAtm = (5E18)/1e12; %Tg
molAtm = 5E18*1E3/30; %g/30 = moles of air
InitAtm13Cdel = -49/1000;
InitAtm13Cratio = (InitAtm13Cdel + 1) * Rstd13C;

InitAtmTMol(1) = S.totalTMol(1) ./ L.k.total; %0.7*1E-6 * molAtm/1E12;
InitAtmTMol(2) = S12C.totalTMol(1) ./ L.k.total; %InitAtmTMol(1)/(1 + InitAtm13Cratio);
InitAtmTMol(3) = S13C.totalTMol(1) ./ (L.k13C.OH + L.k13C.soil + L.k13C.strat); %InitAtmTMol(1)/(1 + (1/InitAtm13Cratio));
InitAtmTMol(4) = SD.totalTMol(1) ./ (L.kD.OH + L.kD.soil + L.kD.strat); 
%check initial atm isotope ratio
InitAtm13Cdel = (((InitAtmTMol(3)./InitAtmTMol(2))./Rstd13C)-1)*1000;

CH4_Mtot = InitAtmTMol(1); 
CH4_delC = InitAtm13Cdel;
CH4_delD = InitAtmTMol(4); 

%run ODE solver for methane
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4]);
% [T,Y] = ode45(@methane1boxderiv_3,[1800 2010],InitAtmTMol,options,S,S13C,S12C,SD,L);

% Data.C13.year = [1700 1900:20:2000];
% Data.C13.del = [-49.0 -49.0 -48.8 -48.6 -48.4 -47.8 -47.1]; 
% Data.ppm.year = [1700 1900:20:2000];
% Data.ppm.ppm = 0.001.*[700 750 900 980 1100 1200 1600 1780 1790]; 

% Data.C13.year = [1800  1900 1950 1955   1960  1965   1970  1975   1980  1985   1990   1995   2000   2005  2010];
% Data.C13.del =  [-49.2 -49  -48  -47.95 -47.9 -47.85 -47.8 -47.75 -47.7 -47.6  -47.5  -47.35 -47.2  -47.2 -47.3];
% Data.ppm.year = [1800 1900 1950 1955 1960 1965 1970 1975 1980 1985 1990 1995 2000 2005 2010];
% Data.ppm.ppm =  [750  850  1050 1150 1200 1300 1400 1450 1580 1660 1720 1740 1770 1760 1780]*0.001;

% figure(1);
% subplot(3,2,1);
%   plot(S.year,S.totalTMol./100,'c',S.year,S13C.totalTMol,'b');
%   title('12C and 13C methane sources');
%   legend('12C/100','13C','Location','NorthWest');
%   legend('boxoff');  
%   %xlim([min(S.year) max(S.year)]);
%   %xlim([1950 2010])
%   ylim([0. .5]);
%   %plot(T,Y(:,2)./100,'c',T,Y(:,3),'b')
%   %title('12C/100 and 13C TMol');
%   %legend('12C/100','13C', 'Location', 'SouthEast');
%   %legend('boxoff')
% subplot(3,2,2);
%   plot(T,Y(:,1).*1E12*1E6 ./molAtm);
%   title('atmospheric methane ppm');
%   %xlim([min(S.year) max(S.year)]);
%   %xlim([1950 2010])
%   %ylim([1.6 2.1]);
%   %set(gca,'YTick',0.5:0.25:2)
%   hold on;
%   plot(Data.ppm.year, Data.ppm.ppm,'o');
%   hold off;
% subplot(3,2,3);
%   plot(S.year,(((S13C.totalTMol./S12C.totalTMol)./Rstd13C)-1)*1000);
%   title('del 13C of methane sources');
%   %xlim([min(S.year) max(S.year)]);
%   %xlim([1950 2010])
%   ylim([-54 -50]);%plot(T,Y(:,3)./Y(:,2))
%   %title('13/12C molar ratio');
% subplot(3,2,4);
%   Ratm=Y(:,3)./Y(:,2); %atmospheric 13C/12C ratio
%   plot(T,((Ratm./Rstd13C) -1)*1000); % plot in del notation
%   title('del 13C atmospheric methane');
%   %xlim([min(S.year) max(S.year)]);
%   %xlim([1950 2010])
%   %ylim([-47 -45.5]);
%   hold on;
%   plot(Data.C13.year, Data.C13.del,'o');
%   hold off;
% subplot(3,2,5);
%   plot(S.year,(((SD.totalTMol./S.totalTMol)./RstdD)-1)*1000);
%   title('del D/H of methane sources');
%   %xlim([1950 2010])
% subplot(3,2,6);
%   Ratm=Y(:,4)./Y(:,2); %atmospheric D/H ratio
%   plot(T,((Ratm./RstdD) -1)*1000); % plot in del notation
%   title('del D atmospheric methane');
%   %xlim([min(S.year) max(S.year)]);  
%   %xlim([1950 2010])
%   ylim([-150 -50]);
%   set(gca,'YTick',-150:20:-50);
%   %figure(10)
%   %compute methane lifetime
%   %plot(T,Y(:,1).*1E12*1E6 ./molAtm);
%   OHSF = Y(:,1)*L.scaling.slope + L.scaling.intercept; %compute OH scaling based on atm methane
%   kTotal = L.k.OH * OHSF + L.k.soil + L.k.strat;
%   %plot(T,kTotal)
%   %title('methane loss - kTotal');