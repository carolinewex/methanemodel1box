% function [source.Tg_micro_CH4]=calc_sources_mircoCH4(I)
ian=1;
source.number=21; %only use number ending in 1
increment=(length(M_tot_int)-1)/(source.number-1);
source.Tg_micro_CH4=1:source.number;
crc_fct=1; %correction factor, how much to change sources for each correction
pct_tol=0.01; %percent tolerance, how much theoret and analyt allowed to vary
for j=1:length(source.Tg_micro_CH4)
while ian>0;
    I.Tg_micro_CH4 = getEm2(tseg,1*source.Tg_micro_CH4, I.tinterp);
    for i = 1:length(I.tinterp)
        I_temp.Tg_micro_CH4 = I.Tg_micro_CH4(i);
    [Out] = calc_emissions(I_temp);
    results.CH4_ppb(:,i) = Out.CH4_ppb;
    int_increment=round((length(results.CH4_ppb)-1)/(source.number-1));
    end
    diff_value =(M_tot_int(1,(increment*j)-(increment-1)) - results.CH4_ppb(1,(int_increment*j)-(int_increment-1)))/results.CH4_ppb(1,(int_increment*j)-(int_increment-1));
    if diff_value <= pct_tol & diff_value >= -pct_tol
        break
    end
    if diff_value > pct_tol
        source.Tg_micro_CH4(1,j)=source.Tg_micro_CH4(1,j)+crc_fct
    end
    if diff_value < (-1*pct_tol)
        source.Tg_micro_CH4(1,j)=source.Tg_micro_CH4(1,j)-crc_fct
    end
end
end
figure (1)
plot(EPICAMethane(:,2),1e-9*mol_atmo*16.04*1e-12*EPICAMethane(:,3),'.k')
hold on
plot(xq,1e-9*mol_atmo*16.04*1e-12*M_tot_int,'-b')
hold on
plot(I.tinterp,1e-9*mol_atmo*16.04*1e-12*results.CH4_ppb,'-g')
hold on
xlabel('Time (years before present)')
xlim([9000 14000])
    set ( gca, 'xdir', 'reverse' )
    xticks(9000:1000:14000)
    ax=gca; ax.XAxis.Exponent = 0;
ylabel('Total Methane (Tg)')
ylim([1200 2100])
hold off