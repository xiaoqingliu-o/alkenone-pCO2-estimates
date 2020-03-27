function output=CO2_Estimate(data)
%
% Pleaes read "read me" to prepare the input data
% Please cite the source publication when using this method:
%
% Zhang et al. (2020).Refining the alkenone-pCO2 method II: towards
% resolving the physiological parameter ?b?.GCA.
%
%INPUTS:
%data = 13 x N vector, including age, 6 parameters and their standard
%devitaions

%OUTPUTS:
%output.Ep = 3 x N vector of inferred Ep, includes 16% level (lower
%1sigma),50% level (median values),84% level (upper 1sigma).
%
%output.b = 3 x N vector of inferred b, includes 16% level (lower
%1sigma),50% level (median values),84% level (upper 1sigma).
%
%output.CO2 = 3 x N vector of inferred CO2, includes 16% level (lower
%1sigma),50% level (median values),84% level (upper 1sigma).
%%
A=-0.4282801;     % A is the slope of the linear regression between Log10 of cell volume and Log10 of growth rate
B=0.574975073;    % B is the intercept of the linear regression between Log10 of cell volume and Log10 of growth rate
f=25;      % enzymatic isotope fractionation associated with intracellular C fixation (per mil)
d=0.7;     % diffusive isotope fractionation of CO2(aq) in seawater (per mil)
     

for j=1:length(data(:,1))
    
    % sample the input parameters within their 2 standard deviations of the
    % mean to obtain 10000000 samples
    T_n=normrnd(data(j,2),data(j,3),[10000000,1]);  
    T_r=T_n(find(T_n>=data(j,2)-2*data(j,3)&T_n<=data(j,2)+2*data(j,3)));
    C37_n=normrnd(data(j,4),data(j,5),[10000000,1]); 
    C37_r=C37_n(find(C37_n>=data(j,4)-2*data(j,5)&C37_n<=data(j,4)+2*data(j,5)));
    Ccarb_n=normrnd(data(j,6),data(j,7),[10000000,1]); 
    Ccarb_r=Ccarb_n(find(Ccarb_n>=data(j,6)-2*data(j,7)&Ccarb_n<=data(j,6)+2*data(j,7)));
    L_n=normrnd(data(j,8),data(j,9),[10000000,1]);
    L_r=L_n(find(L_n>=data(j,8)-2*data(j,9)&L_n<=data(j,8)+2*data(j,9)));    
    S_n=normrnd(data(j,10),data(j,11),[10000000,1]);
    S_r=S_n(find(S_n>=data(j,10)-2*data(j,11)&S_n<=data(j,10)+2*data(j,11)));
    pH_n=normrnd(data(j,12),data(j,13),[10000000,1]);
    pH_r=pH_n(find(pH_n>=data(j,12)-2*data(j,13)&pH_n<=data(j,12)+2*data(j,13)));
    

    P_n=normrnd(5.09,0.08,[10000000,1]);  % cell membrane permeability = (5.09±0.08)*10^(-5) m/s
    P_r=10^(-5).*P_n(find(P_n>=5.09-2*0.08 & P_n<=5.09+2*0.08));
    
    % conduct 10,000 interactions by randomly sampling each input parameter
    % from 10000000 samples
 
    for i=1:10000
        T=T_r(randi(length(T_r),1));            % sea surface temperature (C celcius)
        C37=C37_r(randi(length(C37_r),1));       % alkenone delta 13C value 
        Ccarb=Ccarb_r(randi(length(Ccarb_r),1)); % delta13C of carbonate 
                                                       
        S=S_r(randi(length(S_r),1));      % salinity of seawater
        pH=pH_r(randi(length(pH_r),1));   % pH values of seawater 
        L=L_r(randi(length(L_r),1));      % coccolith length (micrometer)
        P=P_r(randi(length(P_r),1));      % cell membrane permeability to aqueous CO2
        
        E1=11.98-0.12.*T;          % Calcite-CO2(gas) enrichment factor
        E2=-373./(273.15+T)+0.19;   % fractionation factor between gaseous and aqueous CO2
        Cgas=(1000+Ccarb)./(E1/1000+1)-1000; % delta13C of gaseous CO2
        Caq=(1000+Cgas).*(E2./1000+1)-1000;  % delta13C of aqueous CO2
        Corg=(C37+1000).*(4.2/1000+1)-1000;   % delta13C of oragni carbon
        Ep(j,i)=[(Caq+1000)./(Corg+1000)-1]*1000; 
        
        Kw=10.^(-3441./(T+273.15)-2.241+0.09415*sqrt(S));   % ion product of seawater
        k1=8718.*exp(-62.8./(8.3143.*(T+273.15)));       % rate coefficient
        k2=(680.5-4.72.*S)*10^8.*exp(-69400./(8.3143.*(T+273.15))); % rate coefficient
        OH=Kw./(10.^(-pH));     % hydrogen ion concentrations
        k=k1.*OH+k2;            % rate coefficient
        D=5.019*10^(-6).*exp(-19510./(8.3143.*(T+273.15))).*(0.9508-7.389*10^(-4).*T); % diffusivity of aqueous CO2 in seawater
        rk=sqrt(D./k);           % reacto-diffusive length
        r=0.44.*L+0.28;          % cell radius calculated from coccolith mean length
        logV=log10(4/3*pi.*r^3); 
        logu=A.*logV+B;          % Log10 of growth rate calculated using the relationship between cell volume and growth rate
        SE=sqrt(0.09657*(1/89+(logV-1.77276)^2/89/0.0607));
        logu_n=normrnd(logu,SE,[1000,1]);
        logu_r=logu_n(find(logu_n>=logu-2*SE & logu_n<=logu+2*SE));
        Logu= logu_r(randi(length(logu_r),1)); 
        u=10.^(Logu);       
        Q=1.46*10^(-14)*(4/3*pi*r.^3); % cell carbon content derived from cell volume
        b(j,i)=1000*(f-d).*Q.*(u./log(2))/(24*3600)./(4*pi.*(r*10^(-6)).^2).*[r*10^(-6)./(D.*(1+r*10^(-6)./rk))+1./P];

        if f-Ep(j,i)>0
        CO2=b(j,i)/(f-Ep(j,i));
        else CO2=NaN;
        end
        
        K=exp(-60.2409+93.4517*100/(T+273.15)+23.3585*log((T+273.15)/100)+35*(0.023517-0.023656*(T+273.15)/100+0.0047036*((T+273.15)/100)^2));
        pCO2(j,i)=CO2/K; % applying Henry's law
    end
end

output.Ep=prctile(Ep,[16 50 84],2);
output.b=prctile(b,[16 50 84],2);
output.CO2=prctile(pCO2,[16 50 84],2);
