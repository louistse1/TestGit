function [ delu delT cv_stor P_max ] = PREOS( fluidProps, rho_stor, T2, T1 )

%% Peng Robinson for initial cv_stor value
R = 8.314462175;                        % gas constant [=] J/(mol K)
Tc = fluidProps.Tc;
Pc = fluidProps.Pc;
M = fluidProps.M;
w = fluidProps.w;
A3 = fluidProps.A3;  A2 = fluidProps.A2; A1 = fluidProps.A1; A0 = fluidProps.A0;

T_stor_PR = [T2 T1];

for ii = 1:2
    
T = T_stor_PR(ii) + 273.15;

% Reduced variables
Tr = T/Tc;
vm = M/rho_stor;

% Parameters of the EOS for a pure component
kappa = 0.37464 + 1.54226*w - 0.26992*w^2;
alpha = (1 + kappa*(1 - sqrt(Tr)))^2;
a = 0.45724*(R*Tc)^2/Pc*alpha;
b = 0.0778*R*Tc/Pc*1;
ac = 0.45724*(R*Tc)^2/Pc;
bc = 0.0778*R*Tc/Pc;

%Pitzer relation
f0 = 5.92714 - 6.09648/Tr - 1.28862*log(Tr) + 0.169347*Tr^6;
f1 = 15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr^6;
Prsat = exp(f0 + w*f1);     
Psat = Prsat*Pc;            %[kPa] Saturated pressure

Ac = ac*Pc/(R*Tc)^2;
Bc = bc*Pc/(R*Tc);

Asat = a*Psat/(R*T)^2;  
Bsat = b*Psat/R/T;

% Compressibility factor
Zvals_sat = roots([1 -(1-Bsat) (Asat-3*Bsat^2-2*Bsat) -(Asat*Bsat-Bsat^2-Bsat^3)]);
Zvals_c = roots([1 -(1-Bc) (Ac-3*Bc^2-2*Bc) -(Ac*Bc-Bc^2-Bc^3)]);

ZR_sat = [];

for i = 1:3

   if isreal(Zvals_c(i))
   	Zc = Zvals_c(i);
   end
   
   if isreal(Zvals_sat(i))
   	ZR_sat = [ZR_sat Zvals_sat(i)];
   end
   
end

Z_Lsat = min(ZR_sat);
Z_Gsat = max(ZR_sat);

rho_c = Pc/(Zc*R/M*Tc);
rho_Gsat = Psat/(Z_Gsat*R/M*T);
rho_Lsat = Psat/(Z_Lsat*R/M*T);

%% Enthalpy departure function

sq = sqrt(2);

if rho_stor > rho_Gsat && rho_stor < rho_Lsat    %Two-phase
    
    P(ii) = Psat;
    xx = (1/rho_stor - 1/rho_Lsat)/(1/rho_Gsat - 1/rho_Lsat);
    
    A = a*P(ii)/(R*T)^2;
    B = b*P(ii)/(R*T);
    Z = roots([1 -(1-B) (A-3*B^2-2*B) -(A*B-B^2-B^3)]);
   
    dadT = -0.45724*R^2*Tc^2/Pc*kappa*sqrt(alpha/(T*Tc));
    h_Dep_G = 1/(R*Tc)*(R*T*(1-Z_Gsat)-((P(ii)*dadT - R^2*T*Asat)/(2*sq*R*B))*log((Z_Gsat+2.414*B)/(Z_Gsat-0.414*B)));
    h_Dep_L = 1/(R*Tc)*(R*T*(1-Z_Lsat)-((P(ii)*dadT - R^2*T*Asat)/(2*sq*R*B))*log((Z_Lsat+2.414*B)/(Z_Lsat-0.414*B)));
    h_Dep = h_Dep_G*xx + h_Dep_L*(1-xx); 
    
else  %Single-phase
    
    P(ii) = R*T/(vm-b) - ac*alpha/(vm^2 + 2*b*vm - b^2);
    
    A = a*P(ii)/(R*T)^2;
    B = b*P(ii)/(R*T);
    Z = P(ii)*vm/(R*T);
    

    dadT = -0.45724*R^2*Tc^2/Pc*kappa*sqrt(alpha/(T*Tc));
    h_Dep = 1/(R*Tc)*(R*T*(1-Z)-((P(ii)*dadT - R^2*T*A)/(2*sq*R*B))*log((Z+2.414*B)/(Z-0.414*B)));
    
    if P(ii) > Pc && T > Tc     %Supercritical
            xx = 2;
    else
        if rho_stor < rho_Gsat   %Superheated vapor
            xx = 1;
        else
            xx = 0;          %Subcooled liquid
        end
    end
    
end

% Critical enthalpy departure
dadT_c = -0.45724*R^2*Tc^2/Pc*kappa*sqrt(1/(Tc*Tc));
h_Depc = 1/(R*Tc)*(R*Tc*(1-Zc)-((Pc*dadT_c - R^2*Tc*Ac)/(2*sq*R*Bc))*log((Zc+2.414*Bc)/(Zc-0.414*Bc)));

% Ideal gas enthalpy
delh_ideal = A3*(T^4-Tc^4)/4 + A2*(T^3-Tc^3)/3 + A1*(T^2-Tc^2)/2 + A0*(T-Tc);

% Enthalpy and internal energy relative to reference at critical point
h = (R/M*Tc*(h_Depc - h_Dep) + delh_ideal);
u(ii) = h - P(ii)/rho_stor + Pc/rho_c;
end

delu = max(u) - min(u);  %[kJ/kg]
delT = max(T_stor_PR) - min(T_stor_PR);  %[K]
cv_stor = delu/delT*1000;  %[J/kgK]
P_max = P(1)/1000;  %[MPa]