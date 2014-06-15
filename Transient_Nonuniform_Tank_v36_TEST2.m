%%%%%System Model of a Supercritical Thermal Energy Storage System
%%%%%Written by Dr. Adrienne Lavine and Louis Tse
%%%%%Dept. of Mechanical & Aerospace Engineering
% function [ m_stor Ex ] = initscalar( griddata, rundata )
clc; close all;
clearvars
set(0,'DefaultAxesFontSize',14)
%test
%test2
%% Version history
%%v27 additions:
    %Added native MATLAB code for Peng-Robinson Equation of State (no more
    %lookup tables)
%%v26 additions:  
    %Added ability to vary T_stor_initial, in addition to rho_stor
%%v25 additions:  
    %Added exergy analysis section
%%v24 additions:  
    %Added bypass automation so bypass ends exactly at end of discharge (and thus, T_HTF,out = 390C exactly at end of discharge)
    %Added storage term
    
%% Procedure
%(1) Guess a value for the m_stor, the mass of storage fluid.
%(2) Main routine will calculate T_stor(x,t) and T_HTF(x,t).
%(3) T_stor(x,t) at beginning and final time are pinpointed on Peng-Robinson
     %look-up table to determine u_initial and u_final.
%(4) E_stor = m_storage*delu, where E_stor and delu are known.  m_storage is calculated, and
     %iterated upon to match m_stor and m_storage.

%% System specifications


storFluid = 'Naphthalene';
powerBlockconfig = 'config6';
m_stor = 2.4943e7;                   %[kg] Initial guess of m_stor

dischargeTime = 12*3600;        %[s]  Total discharge time 
chargeTime_init = 6*3600;
chargeTime_final = 18*3600;        %[s]  Total charge time 
chargeTime = (chargeTime_final - chargeTime_init)/3600;

T_stor_initial_array = 400:100:500;      %[deg C] (aka T_high)
rho_stor_array = 600;         %[kg/m^3]  Storage fluid loading 
m_dot_charge_array = 100:25:700;
xi_array = 10;


z = 1;                          %Counter for the Results array
z_tot = length(T_stor_initial_array) + length(rho_stor_array) + length(m_dot_charge_array) + length(xi_array);

%% Fluid properties 

%%%Heat transfer fluid properties [1]
rho_HTF = 1060;                         %[kg/m^3]  Density of HTF
cp_HTF = 2446;                          %[J/kgK]  Specific heat of HTF
k_HTF = 0.098;                          %[W/mK] Thermal conductivity of HTF
alpha_HTF = k_HTF/(rho_HTF*cp_HTF);     %[m^2/s]  Thermal diffusivity of HTF
nu_HTF = 9.9e-7;                        %[m^2/s]  Kinematic viscosity of HTF
price_HTF = 3.96;                       %[$/kg]   Price of HTF


%%%Storage fluid properties    
[ Tc Pc w M k_stor price_stor A3 A2 A1 A0 ] = fluidProps(storFluid);
fluidProps.Tc = Tc;
fluidProps.Pc = Pc;
fluidProps.w = w;
fluidProps.M = M;
fluidProps.A3 = A3; fluidProps.A2 = A2; fluidProps.A1 = A1; fluidProps.A0 = A0;

[ m_dot_e, T_powerblock_in, T_powerblock_out, K1, K2, K3, K4, nu ] = powerBlock( powerBlockconfig );
        
%% Main routine

for s4 = 1:length(xi_array)
    xi = xi_array(s4);

 for s3 = 1:length(m_dot_charge_array)
    m_dot_charge = m_dot_charge_array(s3);

  for s2 = 1:length(T_stor_initial_array)
    T_stor_initial = T_stor_initial_array(s2);

    for s = 1:length(rho_stor_array)  %******Put varying parameter here...incremental increase is at the end of this WHILE loop 
      rho_stor = rho_stor_array(s);


    clear bypassZeroes m_dot_b
    bypassZeroes = 3;                           %[h]  Duration of time that the bypass loop is turned off
    m_dot_b_d(1) = 4;                             %Needed to trigger WHILE loop
    r = 1;                                      %Counter for iterations needed to converge to m_stor corresponding to bypass ending exactly at end of discharge
    % T_stor_norm = 11;
    % rr = 1;
    
%     while abs(T_stor_norm) > 10
    while bypassZeroes > 4 || m_dot_b_d(end) > 3   %While loop executes when either condition is TRUE
    
        %The goal of this big IF statement is to iterate to find the value
        %of m_stor for bypass to terminate at exactly end of discharge period.
        
        if r > 1    
            if bypassZeroes > 0                                     %%If bypass is turned off...
                if bypassZeroes < 4                                 %%If bypass is turned off for 2 timesteps, leave m_stor as is
                    m_stor = m_stor;
                elseif bypassZeroes < 6 && bypassZeroes > 3         %%If bypass is turned off for 5 timesteps, increase m_stor 0.04%
                    m_stor = 1.0001*m_stor;
                elseif bypassZeroes < 11 && bypassZeroes > 5        %%If bypass is turned off for 6-10 timesteps, increase m_stor 0.08%
                    m_stor = 1.0007*m_stor;
                elseif bypassZeroes < 26 && bypassZeroes > 10       %%If bypass is turned off for 11-25 timesteps, increase m_stor 0.1%
                    m_stor = 1.0025*m_stor;
                elseif bypassZeroes < 51 && bypassZeroes > 25       %%If bypass is turned off for 26-50 timesteps, increase m_stor 0.3%
                    m_stor = 1.0035*m_stor;
                elseif bypassZeroes < 101 && bypassZeroes > 50      %%If bypass is turned off for 51-100 timesteps, increase m_stor 0.5%
                    m_stor = 1.01*m_stor;
                elseif bypassZeroes < 200 && bypassZeroes > 100     %%If bypass is turned off for 101-200 timesteps, increase m_stor 1%
                    m_stor = 1.0155*m_stor;
                elseif bypassZeroes < 300 && bypassZeroes > 201     %%If bypass is turned off for 101-200 timesteps, increase m_stor 1%
                    m_stor = 1.0235*m_stor;                    
                else                                                %%If bypass is turned off for more than 200 time steps, increase m_stor by 10%
                    m_stor = 1.035*m_stor;
                end
            else                                                    %%If bypass is turned on...just decrease m_stor so much that it goes into a no-bypass mode
                if m_dot_b_d(end) < 3    
                    m_stor = m_stor;
                elseif m_dot_b_d(end) < 10 && m_dot_b_d(end) >= 3
                    m_stor = 0.9985*m_stor;    
                elseif m_dot_b_d(end) < 20 && m_dot_b_d(end) >= 10
                    m_stor = 0.992*m_stor;
                elseif m_dot_b_d(end) < 50 && m_dot_b_d(end) >= 20
                    m_stor = 0.99*m_stor;
                elseif m_dot_b_d(end) < 100 && m_dot_b_d(end) >= 50
                    m_stor = 0.985*m_stor;  
                elseif m_dot_b_d(end) < 150 && m_dot_b_d(end) >= 100
                    m_stor = 0.95*m_stor; 
                elseif m_dot_b_d(end) < 250 && m_dot_b_d(end) >= 150
                    m_stor = 0.9*m_stor; 
                else
                    m_stor = 0.8*m_stor;
                end               
            end
        end

clear cv_stor

%% Peng Robinson for initial cv_stor value
T_stor_initial = T_stor_initial_array(s2);
rho_stor = rho_stor_array(s);
[ delu delT cv_stor_initial P_max ] = PREOS( fluidProps, rho_stor, T_stor_initial, 0.999*T_stor_initial );

%% Tube properties

    %%%Tube material properties
    F_tu = 505;         %[MPa]  Ultimate tensile strength of SS316L
    n = 4;              %Safety factor
    d = 0.6;            %Derating factor for SS316L operating up to 500 deg C
    sigma_SS = F_tu*d/n;   %[MPa]  Allowable stress
    price_SS = 1.40;    %[$/kg]  Price of tube material
    cp_SS = 500;        %[J/kgK] Specific heat of tube material
    rho_SS = 7990;      %[kg/m^3] Density of tube material
    k_SS = 16.3;        %[W/mK] Thermal conductivity of tube material
    
    %%%Tube dimensions
    r1 = 0.025;        %[m]  Radius of a single tube, with zero wall thickness  
    L = 30;             %[m]  Length of a single tube            
    h_i = 18;            %[W/m^2K]  Overall heat transfer coefficient from storage fluid-to-steel
    h_o = 500;            %[W/m^2K]  Overall heat transfer coefficient from steel-to-HTF
    U_tot = (1/h_i + 1/h_o)^-1;
    t_wall = P_max*r1/sigma_SS;  %[m] Required wall thickness    
    r2 = r1 + t_wall;      %[m]  Outer radius of tube
    r3 = r2 + 0.001;        %[m]  Annulus of HTF surrounding a single tube, which implies Gap = (r_o - r1)    
        
    %%%Tube bundle geometry
    S_T = 2.25*r1;            %[m] Transverse pitch
    S_L = 2.25*r1;            %[m] Longitudinal pitch
    S_D = (S_L^2 + (S_T/2)^2)^0.5;  %[m] Diagonal pitch
    
    if S_T < 2*r1 || S_L < 2*r1
        break
        print('Tube geometry interference error')
    end
    
%% Set grid and time steps

        timeNodes_d = 9000;              
        dt_d = dischargeTime/timeNodes_d;                %[s]
        time_d = 0:dt_d/3600:dischargeTime/3600;          %[h]
        t_final_d = max(size(time_d));      %index of final time node

        lengthNodes_d = 9000;
        dx_d = L/lengthNodes_d;
        tankLength_d = 0:dx_d:L;
        x_final_d = max(size(tankLength_d));       %index of final length node

        %%%Matrix pre-allocation%%%
        T_stor_d = zeros(lengthNodes_d+1, timeNodes_d+1);  %Pre-allocate matrix to decrease computing time
        T_SS_d = zeros(lengthNodes_d+1, timeNodes_d+1);   %Pre-allocate matrix to decrease computing tim
        T_HTF_d = zeros(lengthNodes_d+1, timeNodes_d+1);   %Pre-allocate matrix to decrease computing time
        
        T_stor_avg_d = zeros(timeNodes_d+1, 1);
        T_SS_avg_d = zeros(timeNodes_d+1, 1);
        T_HTF_avg_d = zeros(timeNodes_d+1, 1);
        
        T_ei_d = zeros(timeNodes_d+1, 1);
        T_eo_d = zeros(timeNodes_d+1, 1);
        m_dot_b_d = zeros(timeNodes_d+1, 1);
        m_dot_s_d = zeros(timeNodes_d+1, 1);
        
%% Set initial and inlet conditions%%%

        T_ei_d(1) = T_powerblock_in;             %[deg C] (design value of turbine)
        T_eo_d(1) = K1*T_ei_d(1) + K2;             %[deg C] (aka T_low)

        % if rr == 1
            
        %     T_stor_d(:,1) = T_stor_initial;
        %     T_stor_avg_d(1) = T_stor_initial;
        %     T_SS_d(:,1) = T_stor_initial;
        %     T_SS_avg_d(1) = T_stor_initial;
        %     T_HTF_d(1,1) = T_eo_d(1);
        % else
        %     T_stor_d(:,1) = T_stor_d_initial(:,1);
        %     T_stor_avg_d(1) = T_stor_avg_d_initial(1);
        %     T_SS_d(:,1) = T_SS_d_initial(:,1);
        %     T_SS_avg_d(1) = T_SS_avg_d_initial(1);
        %     T_HTF_d(1,1) = T_eo_d(1);
        % end
        
        m_dot_s_d(1) = 261.8341;                 %[kg/s]  (total guess, but is needed to start the routine)  
        m_dot_b_d(1) = m_dot_e - m_dot_s_d(1);  %[kg/s]  Mass flow rate through the bypass (to cool down the hot HTF exiting the tank)
        
%% Calculated values

        %%%Calculated values%%%
        N = m_stor/(pi*r1^2*rho_stor*L);    %Number of tubes
        P_h1 = 2*pi*r1*N;                  %[m]  Total perimeter of storage fluid in contact with steel
        P_h2 = 2*pi*r2*N;            %[m]  Total perimeter of steel in contact with HTF
        
        Ac_stor = pi*r1^2*N;                   %[m^2]    Total cross-sectional area of storage fluid
        Ac_SS = pi*(r2^2 - r1^2)*N;             %[m^2]    Total cross-sectional area of steel
        Ac_HTF = pi*(r3^2 - r2^2)*N;          %[m^2]  Total cross-sectional area of HTF
        
        m_prime_stor = rho_stor*Ac_stor;        %[kg/m]    Mass of storage fluid per length of all tubes
        m_prime_SS = rho_SS*Ac_SS;
        
        V_stor = pi*(r1)^2*L*N;                %[kg/m^3]  Volume of storage fluid
        u_m = m_dot_e/(rho_HTF*Ac_HTF);         %[m/s]
      
        R1_stor = h_i*P_h1*dt_d/(m_prime_stor);
        R1_SS = h_o*P_h2*dt_d/(m_prime_SS*cp_SS);
        R2_SS = h_i*P_h1*dt_d/(m_prime_SS*cp_SS);
        R3_SS = k_SS*Ac_SS/(m_prime_SS*cp_SS*dx_d^2);
        R1_HTF = h_o*P_h2*dx_d/cp_HTF;

             
        
%% DISCHARGING PERIOD

tic;
tstart = tic;

for t = 1:timeNodes_d
    
    clear delu delT delu_delT cv_stor_tracking cv_stor_error
    q = 1;                          %Counter for iteratiosn needed to converge to correct cv_stor
    cv_stor(t) = 2233;                 %[J/kgK]  Initial guess of specific heat of storage fluid
    cv_stor_error(1) = 1;           %Needed to trigger WHILE loop
        
    while abs(cv_stor_error) > 0.01
        
        
    %%%MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        R2_HTF = R1_HTF/m_dot_s_d(t);
        R2_stor = R1_stor/cv_stor(t);
        
       %%%FIRST TIME-STEP (To simplify code...the first time step has isothermal BC's for axial conduction, and no incoming HTF)
       if t == 1  

           for x = 1:lengthNodes_d     
                T_HTF_d(x+1,t) = T_HTF_d(x,t) + R2_HTF*(T_SS_d(x,t) - T_HTF_d(x,t));
           end
           T_SS_d(:,t+1) = T_SS_d(:,t) - R1_SS*(T_SS_d(:,t) - T_HTF_d(:,t)) + R2_SS*(T_stor_d(:,t) - T_SS_d(:,t));
           T_stor_d(:,t+1) = T_stor_d(:,t) - R2_stor*(T_stor_d(:,t) - T_SS_d(:,t));

       %%%REMAINING TIME-STEPS    
       elseif t > 1 

           T_HTF_d(1,t) = T_eo_d(t-1);     %Incoming HTF from the evaporator

           for x = 1:lengthNodes_d     
                T_HTF_d(x+1,t) = T_HTF_d(x,t) + R2_HTF*(T_SS_d(x,t) - T_HTF_d(x,t));
           end
           T_SS_d(:,t+1) = T_SS_d(:,t) - R1_SS*(T_SS_d(:,t) - T_HTF_d(:,t)) + R2_SS*(T_stor_d(:,t) - T_SS_d(:,t));
           T_stor_d(:,t+1) = T_stor_d(:,t) - R2_stor*(T_stor_d(:,t) - T_SS_d(:,t));

       end

        %Average temperatures
        T_stor_avg_d(t+1) = mean(T_stor_d(:,t+1));  %T_stor_avg is used to calculate u_initial and u_final
        T_SS_avg_d(t+1) = mean(T_SS_d(:,t+1));  %T_stor_avg is used to calculate u_initial and u_final
        T_HTF_avg_d(t) = mean(T_HTF_d(:,t));
    
    %%%BYPASS LOGIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %Bypass loop is on.  This means T_ei is set, and m_dot_b is solved for.
        if T_HTF_d(x_final_d,t) > T_powerblock_in

            T_ei_d(t) = T_powerblock_in; 
            T_eo_d(t) = 0.433*T_ei_d(t) + 120.13;    %Kolb function [1]

            m_dot_b_d(t+1) = m_dot_e*(T_ei_d(t)-T_HTF_d(x_final_d,t))/(T_eo_d(t)-T_HTF_d(x_final_d,t));   %Conservation of enthalpy
            m_dot_s_d(t+1) = m_dot_e - m_dot_b_d(t);    %Conservation of mass

        %Bypass loop is off.  This means m_dot_b is set to zero, and T_ei is equal to the HTF temp. exiting the tank.
        else

            T_ei_d(t) = T_HTF_d(x_final_d,t);
            T_eo_d(t) = K1*T_ei_d(t) + K2;

            m_dot_b_d(t+1) = 0;
            m_dot_s_d(t+1) = m_dot_e - m_dot_b_d(t);

        end  %End of bypass loop logic

        
    %%%ITERATING ON CV,STOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %% Peng-Robinson

            T2 = T_stor_avg_d(t); 
            T1 = T_stor_avg_d(t+1);
            T_stor_diff = T_stor_avg_d(t) - T_stor_avg_d(t+1);
            
            if T_stor_diff < 2e-3
                T2 = T_stor_avg_d(t);
                T1 = 0.999*T_stor_avg_d(t+1);
            end
                
            [ delu delT cv_stor_new P_max ] = PREOS( fluidProps, rho_stor, T2, T1 );    

        %%%Check error on cv_stor and iterate if needed            
        cv_stor_old = cv_stor(t);
        cv_stor_error = (cv_stor_new - cv_stor_old)/cv_stor_new;
        cv_stor(t) = (cv_stor_old + cv_stor_new)/2;                              %q is the index counter for the number of iterations it takes for one case to converge (q resets to 1 after one case has converged)
        
        q = q+1;
        if mod(q,100) == 0
            disp('cv_stor iteration not converging')
            break
        end
        
    end   %End of looping for correct cv_stor (q is reset to 1)

 end %End of marching forward in time for entire discharge period

%% Power block

    Q_dot_turbine = K3*T_ei_d + K4;      %[MW], Kolb function [1]
    Q_total = sum(Q_dot_turbine*dt_d/3600);   %[MWh]
 
 
%% Clean-up

    %These are to copy an extra step in time for each vector to make it the same length as the Time vector.
    T_HTF_d(:,timeNodes_d+1) = T_HTF_d(:,timeNodes_d);                    
    T_ei_d(t_final_d) = T_ei_d(t_final_d-1);
    T_eo_d(t_final_d) = T_eo_d(t_final_d-1);
    m_dot_b_d(timeNodes_d+1) = m_dot_b_d(timeNodes_d);              
    Q_dot_turbine(t_final_d) = Q_dot_turbine(t_final_d-1);
    cv_stor(t_final_d) = cv_stor(t_final_d-1);
    cv_stor = cv_stor(:);
    
    %This is to copy the last row in x (which is all zeroes) from the next-to-last row
    T_stor_d(lengthNodes_d+1,:) = T_stor_d(lengthNodes_d,:);

    %These are to copy an extra step in x to make the vector the same length as the Length vector.
    T_HTF_d(x_final_d,timeNodes_d+1) = T_HTF_d(x_final_d-1,timeNodes_d+1);  
    
    %This is to remove the last time step
    T_HTF_avg_d = T_HTF_avg_d(1:t_final_d-1);
    

%% Energy and exergy analysis

    %Energy and exergy analyses
    T0 = 25 + 273.15;
    dE_out = m_dot_e*cv_stor(:).*(T_ei_d(:)-T_eo_d(:))*dt_d;
    dEx_out = m_dot_e*cv_stor(:).*(T_ei_d(:)-T_eo_d(:) - T0.*log((T_ei_d(:)+273.15)./(T_eo_d(:)+273.15)))*dt_d;
    dQ_loss = (0.00017*T_stor_avg_d(:)+0.012)*dt_d;       %[kJ/m^2];

    E_out_tot = sum(dE_out(:))/(1e6*3600);        %[MWh]
    Ex_out_tot = sum(dEx_out(:))/(1e6*3600);    %[MWh]   
    
    
    
    %% Peng Robinson for final energy calculations
    
        T2 = T_stor_avg_d(1);
        T1 = T_stor_avg_d(end);
        [ delu_final delT cv_stor_final ~ ] = PREOS( fluidProps, rho_stor, T2, T1 ); 

        %%%Temperatures
        delT_stor = T_stor_avg_d(1) - T_stor_avg_d(end);
        delT_HTF = T_HTF_avg_d(1) - T_HTF_avg_d(end);

    %% Cost Analysis

        %%%Steel tube parameters
        V_SS = pi*(r2^2 - r1^2)*L*N;   %[m^3]  Total volume of tube material
        m_SS = rho_SS*V_SS;              %[kg]  Total mass of tube material

        %%%HTF parameters
        V_HTF = pi*(r3^2 - r2^2)*L*N;     %[m^3] Total volume of HTF
        m_HTF = rho_HTF*V_HTF;              %[kg]  Total mass of HTF

        %%%Calculating costs
        cost_SS = m_SS*price_SS/(Q_total*1e3/.37);          %[$/kWh]   Cost of steel tube material
        cost_HTF = m_HTF*price_HTF/(Q_total*1e3/0.37);      %[$/kWh]   Cost of HTF
        cost_stor = m_stor*price_stor/(Q_total*1e3/.37);    %[$/kWh]  Cost of storage fluid
        cost_total = cost_HTF + cost_SS + cost_stor;        %[$/kWh]  Total cost

        %%%Energy storage breakdown
        E_SS = m_SS*cp_SS*max(delT_stor)/(1e3*3600);         %[kWh]  Total energy stored in steel
        E_HTF = m_HTF*cp_HTF*max(delT_HTF)/(1e3*3600);  %[kWh]  Total energy stored in HTF
        E_stor = m_stor*delu_final/3600;                %[kWh]  Total energy stored in storage fluid
        E_total = (E_SS + E_HTF + E_stor)/1000;         %[MWh]  Total energy stored

        E_HTF_frac = E_HTF/(E_HTF + E_SS + E_stor);
        E_SS_frac = E_SS/(E_HTF + E_SS + E_stor);
        E_stor_frac = E_stor/(E_HTF + E_SS + E_stor);
            
        
        mc_stor = m_prime_stor*cv_stor_initial;
        mc_SS = m_prime_SS*cp_SS;
        

%%     
    
    
    
    
    
    
    
    
    
 
 
 
%% CHARGING PERIOD
 
clear cv_stor
T_stor_initial_charge = 311;      %[deg C] 
E_discharge = 1621*(1e6*3600);         %[J] converted from 1621 MWh


 %% Set grid and time steps

        timeNodes_c = lengthNodes_d*4;              
        dt_c = (chargeTime_final - chargeTime_init)/timeNodes_c;                %[s]
        time_c = chargeTime_init/3600:dt_c/3600:chargeTime_final/3600;          %[h]
        t_final_c = max(size(time_c));      %index of final time node

        lengthNodes_c = lengthNodes_d;
        dx_c = L/lengthNodes_c;
        tankLength_c = 0:dx_c:L;
        x_final_c = max(size(tankLength_c));       %index of final length node

        %%%Matrix pre-allocation%%%
        T_stor_c = zeros(lengthNodes_c+1, timeNodes_c+1);  %Pre-allocate matrix to decrease computing time
        T_SS_c = zeros(lengthNodes_c+1, timeNodes_c+1);   %Pre-allocate matrix to decrease computing time
        T_HTF_c = zeros(lengthNodes_c+1, timeNodes_c+1);   %Pre-allocate matrix to decrease computing time
        
        T_stor_avg_c = zeros(timeNodes_c+1, 1);
        T_SS_avg_c = zeros(timeNodes_c+1, 1);
        T_HTF_avg_c = zeros(timeNodes_c+1, 1);
        
        
        T_in_c = zeros(timeNodes_c+1, 1);
        T_out_c = zeros(timeNodes_c+1, 1);
        
%% Set initial and inlet conditions%%%

        T_stor_c(:,1) = T_stor_d(:,end);
        T_stor_avg_c(1) = T_stor_avg_d(end);
        T_HTF_c(:,1) = T_HTF_d(:,end);
        T_SS_c(:,1) = T_SS_d(:,end);
        T_SS_avg_c(1) = T_SS_avg_c(end);
        
        m_dot_s_c = m_dot_charge;  
        m_dot_b_c = 0;
 
%% Solar field specs and determining A_solar

        for t = 1:timeNodes_c
            Q_solar(t) = max(917*sin(0.22*time_c(t)-7.51),0);         %[W/m^2] time in Hours, Fit from Daggett, CA weather data
            dQ(t) = Q_solar(t)*dt_c;                %[J/m^2]
        end
        Q_solar_tot = sum(dQ(:));       %[J/m^2];
        A_solar = E_discharge/Q_solar_tot;      %[m^2] Area needed to have energy input = 1621 MWh

%% Calculated values

        %%%Calculated values%%%
        N = m_stor/(pi*r1^2*rho_stor*L);    %Number of tubes
        P_h1 = 2*pi*r1*N;                  %[m]  Total perimeter of storage fluid in contact with steel
        P_h2 = 2*pi*r2*N;            %[m]  Total perimeter of steel in contact with HTF
        
        Ac_stor = pi*r1^2*N;                   %[m^2]    Total cross-sectional area of storage fluid
        Ac_SS = pi*(r2^2 - r1^2)*N;             %[m^2]    Total cross-sectional area of steel
        Ac_HTF = pi*(r3^2 - r2^2)*N;          %[m^2]  Total cross-sectional area of HTF
        
        m_prime_stor = rho_stor*Ac_stor;        %[kg/m]    Mass of storage fluid per length of all tubes
        m_prime_SS = rho_SS*Ac_SS;
        
        V_stor = pi*(r1)^2*L*N;                %[kg/m^3]  Volume of storage fluid
        u_m = m_dot_e/(rho_HTF*Ac_HTF);         %[m/s]
      
%         R1 = U2*P_h2*dx/cp_HTF;                          %Define constants that can be overwritten during loops to decrease computing time
%         R2 = U1*P_h1*dt/(m_prime_stor);          %Define constants that can be overwritten during loops to decrease computing time
%         R5 = dx/(dt*u_m);

        R1_stor = h_i*P_h1*dt_c/(m_prime_stor);
        R1_SS = h_o*P_h2*dt_c/(m_prime_SS*cp_SS);
        R2_SS = h_i*P_h1*dt_c/(m_prime_SS*cp_SS);
        R3_SS = k_SS*Ac_SS/(m_prime_SS*cp_SS*dx_c^2);
        R1_HTF = h_o*P_h2*dx_c/cp_HTF;   
        
%% Main routine


for t = 1:timeNodes_c
    
    clear delu delT delu_delT cv_stor_tracking cv_stor_error
    q = 1;                          %Counter for iteratiosn needed to converge to correct cv_stor
    cv_stor(t) = 2233;                 %[J/kgK]  Initial guess of specific heat of storage fluid
    cv_stor_error(1) = 1;           %Needed to trigger WHILE loop
        
    while abs(cv_stor_error) > 0.01
        
        
    %%%MAIN ROUTINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        R2_HTF = R1_HTF/m_dot_s_c;
        R2_stor = R1_stor/cv_stor(t);
                
       %%%FIRST TIME-STEP (To simplify code...the first time step does not have axial conduction, and no incoming HTF)
       if t == 1  

           for x = 1:lengthNodes_c     
                T_HTF_c(x+1,t) = T_HTF_c(x,t) + R2_HTF*(T_SS_c(x,t) - T_HTF_c(x,t));
           end
           T_SS_c(:,t+1) = T_SS_c(:,t) - R1_SS*(T_SS_c(:,t) - T_HTF_c(:,t)) + R2_SS*(T_stor_c(:,t) - T_SS_c(:,t));
           T_stor_c(:,t+1) = T_stor_c(:,t) - R2_stor*(T_stor_c(:,t) - T_SS_c(:,t));
           
       %%%REMAINING TIME-STEPS    
       elseif t > 1 

           T_HTF_c(1,t) = T_in_c(t-1);     %Incoming HTF from the solar field

           for x = 1:lengthNodes_c     
                T_HTF_c(x+1,t) = T_HTF_c(x,t) + R2_HTF*(T_SS_c(x,t) - T_HTF_c(x,t));
           end
           T_SS_c(:,t+1) = T_SS_c(:,t) - R1_SS*(T_SS_c(:,t) - T_HTF_c(:,t)) + R2_SS*(T_stor_c(:,t) - T_SS_c(:,t));
           T_stor_c(:,t+1) = T_stor_c(:,t) - R2_stor*(T_stor_c(:,t) - T_SS_c(:,t));

       end

        %Average temperatures
        T_stor_avg_c(t+1) = mean(T_stor_c(:,t+1));  %T_stor_avg is used to calculate u_initial and u_final
        T_SS_avg_c(t+1) = mean(T_SS_c(:,t+1));
        T_HTF_avg_c(t) = mean(T_HTF_c(:,t));
    
        %Energy and exergy analyses
        dE_in(t) = Q_solar(t)*A_solar*dt_c;
        dEx_in(t) = (T_HTF_c(1,t)-T_HTF_c(end,t))*dt_c;
        dEx_HT(t) = log((T_HTF_c(1,t)+273.15)/(T_HTF_c(end,t)+273.15))*dt_c;
        dQ_loss(t) = (0.00017*T_stor_avg_c(t)+0.012)*dt_c;       %[kJ/m^2];

        T_in_c(t) = Q_solar(t)*A_solar/(m_dot_s_c*cp_HTF) + T_HTF_c(end,t);
        
      
    %%%ITERATING ON CV,STOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %% Peng-Robinson

            T2 = T_stor_avg_c(t); 
            T1 = T_stor_avg_c(t+1);
            T_stor_diff = T_stor_avg_c(t) - T_stor_avg_c(t+1);
            
            if T_stor_diff < 2e-3
                T2 = T_stor_avg_c(t);
                T1 = 0.999*T_stor_avg_c(t+1);
            end
                
            [ delu delT cv_stor_new P_max ] = PREOS( fluidProps, rho_stor, T2, T1 );    

        %%%Check error on cv_stor and iterate if needed            
        cv_stor_old = cv_stor(t);
        cv_stor_error = (cv_stor_new - cv_stor_old)/cv_stor_new;
        cv_stor(t) = (cv_stor_old + cv_stor_new)/2;                              %q is the index counter for the number of iterations it takes for one case to converge (q resets to 1 after one case has converged)
        
        q = q+1;
        if mod(q,100) == 0
            disp('cv_stor iteration not converging')
            break
        end
        
    end   %End of looping for correct cv_stor (q is reset to 1)

 end %End of marching forward in time for entire discharge period

    
%% Pressure drop
    
    xi = xi_array(s4);
    fricFac = 12;
    N_L = 60;           %[] Number of tube rows
    N_C = N/N_L;        %[] Number of tubes in each row
    Ac_HTF_2 = N_L*(2*r1) + (N_L-1)*S_T;
    u_m = m_dot_charge/(rho_HTF*Ac_HTF_2);         %[m/s]

    if S_D < (S_T + 2*r1)
        u_max = S_T/(2*(S_D - 2*r1))*u_m;
    else
        u_max = S_T/(S_T - 2*r1)*u_m;
    end

    Re_max = u_max*2*r1/nu_HTF;
    deltaP = N_C*xi*(rho_HTF*u_max^2/2)*fricFac/1e3; %[kPa]
    nu_pump = 0.6;
    E_pump_c = m_dot_charge*(chargeTime_final-chargeTime_init)*deltaP/(rho_HTF*nu_pump)/(1e3*3600);  %[MWh] 
    E_pump_d = m_dot_s_d(:)*deltaP/(rho_HTF*nu_pump)*dt_d/(1e3*3600);  %[MWh] 
    E_pump_d = sum(E_pump_d(:));
    E_pump = E_pump_c + E_pump_d;
 


%% Clean-up 
    
    T_HTF_c(:,timeNodes_c+1) = T_HTF_c(:,timeNodes_c);                    
    T_in_c(timeNodes_c+1) = T_in_c(timeNodes_c);           
    T_HTF_avg_c(timeNodes_c+1) = T_HTF_avg_c(timeNodes_c);
    Q_solar(timeNodes_c+1) = Q_solar(timeNodes_c);
    cv_stor(timeNodes_c+1) = cv_stor(timeNodes_c);

    %This is to copy the last row in x (which is all zeroes) from the next-to-last row
    T_stor_c(lengthNodes_c+1,:) = T_stor_c(lengthNodes_c,:);
    time_c = time_c(:);
    
   %% Cost Analysis

    results_T_HTF_avg(:,s3) = T_HTF_avg_c(:);
    results_T_stor_avg(:,s3) = T_stor_avg_c(:);
    results_T_HTF_in(:,s3) = T_in_c(:); 
    results_T_HTF_out(:,s3) = T_HTF_c(end,:);

    %%Energy
    E_in_total_c = sum(dE_in(:))/(1e6*3600);
    firstLawEffic_charge = (E_in_total_c - E_pump)/E_in_total_c;

    %%Exergy
    T0 = 25+273.15;
    T_out(:,1) = T_HTF_c(end,:);
    dEx_in2 = m_dot_charge*cv_stor(:).*(T_in_c(:)-T_out_c(:) - T0.*log((T_in_c(:)+273.15)./(T_out_c(:)+273.15)))*dt_c;
%     dQ_loss = (0.00017*T_stor_avg_c(:)+0.012)*dt;       %[kJ/m^2];

    Ex_in_tot2 = sum(dEx_in2(:))/(1e6*3600);    %[MWh]
    
    dEx_HT_results(s3,:) = dEx_HT;
    Ex_HT_results(s3,:) = m_dot_charge*cp_HTF*T0*dEx_HT_results(s3,:)/(1e6*3600); 

    
    Ex_HT_total = m_dot_charge*cp_HTF*T0*sum(dEx_HT_results(s3,:))/(1e6*3600);               %[MWh]
    Ex_dP_total = 1.25*E_pump;                               %[MWh]
    Ex_in_total = m_dot_charge*cp_HTF*sum(dEx_in(:))/(1e6*3600);
    Ex_stored_total = Ex_in_total - Ex_HT_total;       %[MWh]  Total supplied exergy by the HTF
    Ex_net = Ex_out_tot - Ex_dP_total;     %[MWh] constant for now...
    %Energy and exergy analyses
    
    exEff = Ex_net/Ex_in_tot2;
    
        
    % sumdiff = sum(dEx_in);
    % sumlog = sum(dEx_HT);
    % E_in = m_dot_s*cp_HTF*(sumdiff - T0*sumlog);

    secondLawEffic_charge = Ex_stored_total/Ex_in_total;
    secondLawEffic_total_dP = Ex_net/Ex_stored_total;
    secondLawEffic_total = Ex_out_tot/Ex_stored_total;
    



    
%% Continuity error

    % T_stor_diff = T_stor_d(:,1) - T_stor_c(:,end);
    % T_SS_diff = T_stor_d(:,1) - T_SS_c(:,end);
    % T_HTF_diff = T_stor_d(:,1) - T_HTF_c(:,end);
    % T_stor_norm = sum(T_stor_diff(:));
    % T_SS_norm = sum(T_SS_diff(:));
    % T_HTF_norm = sum(T_HTF_diff(:));
    
    
    % T_stor_avg_diff = T_stor_avg_d(1) - T_stor_avg_c(end);
    % T_SS_avg_diff = T_SS_avg_d(1) - T_SS_avg_c(end);
    % T_HTF_avg_diff = T_stor_avg_d(1) - T_HTF_avg_c(end);

    % %%Save temperatures for new initial conditions
    % T_stor_d_initial(:,1) = T_stor_c(:,end);
    % T_SS_d_initial(:,1) = T_SS_c(:,end);
    % T_HTF_d_initial(:,1) = T_HTF_c(:,end);
    % T_stor_avg_d_initial = T_stor_avg_c(end);
    % T_SS_avg_d_initial = T_SS_avg_c(end);
    
%     checkStatus_cont1 = sprintf('ITERATION %d...T_stor norm = %0.2f',rr,T_stor_norm);
%     checkStatus_cont2 = sprintf('...1)  T_SS norm = %0.2f',T_SS_norm);
%     checkStatus_cont3 = sprintf('...1)  T_HTF norm = %0.2f',T_HTF_norm);
%     disp(checkStatus_cont1)
%     disp(checkStatus_cont2)
%     disp(checkStatus_cont3)
    % rr = rr + 1;
    % end
    %% Automating bypass loop to end exactly at end of discharge period    
    
    %Determine how close bypass is to ending exactly at end of discharge
    bypassZeroes = length(m_dot_b_d) - nnz(m_dot_b_d);                  %Number of zero elements in the bypass mass flow rate array
    bypassTime_off = bypassZeroes*dt_d/3600;                          %[h] Duration of time when bypass is turned off
    bypassTime_on = dischargeTime/3600 - bypassTime_off;            %[h] Duration of time when bypass is turned on
       
    %Print out statuses for user 
    telapsed = toc(tstart);
    
    checkStatus_iter = sprintf('ITERATION %d (%0.1fs)...Case: rho = %d, T_stor = %d',r,telapsed,rho_stor,T_stor_initial);
    checkStatus_mstor = sprintf('...1) Storage fluid mass = %0.4E',m_stor);
    checkStatus_mdotb = sprintf('...2) Final bypass flow rate = %0.2f',m_dot_b_d(end));
    checkStatus_bypassZeroes = sprintf('...3) Bypass zeroes = %d',bypassZeroes);
    
    disp(checkStatus_iter)
    disp(checkStatus_mstor)
    disp(checkStatus_mdotb)
    disp(checkStatus_bypassZeroes)
      
    r = r + 1;
    
    end  %End of automating bypass loop to end exactly at discharge

resultsNumerical_firstLaw(z,:) = [A_solar, E_in_total_c, E_pump, firstLawEffic_charge];
resultsNumerical_secondLaw(z,:) = [m_dot_charge, Ex_HT_total, Ex_dP_total, Ex_in_total, Ex_stored_total, Ex_net...
    secondLawEffic_charge, secondLawEffic_total, secondLawEffic_total_dP, exEff]; 
resultsNumerical_dP(z,:) = [xi, fricFac, N_L, nu_pump];
resultsNumerical_Energy(z,:) = [E_out_tot, E_stor, Ex_out_tot, E_pump, mc_stor, mc_SS];
resultsNumerical_Settings(z,:) = [dischargeTime, N, h_i, h_o, U_tot, L, timeNodes_c, lengthNodes_c, timeNodes_d, lengthNodes_d];
resultsNumerical_Cost(z,:) = [rho_stor, T_stor_initial, t_wall, r2, V_stor, m_stor, E_stor/1000, V_SS, m_SS, E_SS/1000,...
    V_HTF, m_HTF, E_HTF/1000, cost_HTF, cost_SS, cost_stor, cost_total];


    telapsed = toc(tstart);
    checkStatus_runTime = sprintf('ITERATION %d/%d (%0.1fs)',z,z_tot, telapsed);
    disp(checkStatus_runTime)


    z = z+1;
    
 end  %End of looping through different cases (s)
end %End of looping through different cases (s2)
end %s3
end %s4
%% Results array

    header_firstLaw = {'Solar area (m^2)','Energy input (MWh)','Pump work (MWh)','Energy efficiency'};
    header_secondLaw = {'Mass flow rate (kg/s)', 'Exergy loss: HT (MWh)', 'Exergy loss: dP (MWh)', 'Total exergy input (MWh)',...
        'Total exergy stored (MWh)', 'Exergy recovered (MWh)','Exergetic efficiency charging', 'Exergetic effic (no dP)', 'Exergetic effic (dP)', 'exEff'};
    header_dP = {'xi', 'Friction factor', 'Number of tube rows', 'Pump efficiency'};
    header_Energy = {'Energy output (MWh)', 'm*delu (MWh)', 'Exergy output (MWh)', 'Energy pumping (MWh)', 'Heat capacity of stor fluid (J/K)', 'Heat capacity of steel (J/K)'};
    header_Settings = {'Discharge (s)', 'Number of tubes', 'U1 (W/m^2-K)','U2 (W/m^2-K)','U_tot (W/m^2-K)','L (m)', 't nodes charge', 'x nodes charge', 't nodes discharge', 'x nodes discharge'};
    header_Cost = {'Storage fluid loading','T_stor_initial','Tube thickness (m)','Outer tube radius (m)','Storage fluid volume (m^3)','Storage fluid mass (kg)','Energy stored in stor fluid (MWh)','Steel volume (m^3)','Steel mass (kg)','Energy stored in steel (MWh)','HTF volume (m^3)','HTF mass (kg)','Energy stored in HTF (MWh)','HTF cost ($/kWh)','Steel cost ($/kWh)','Storage fluid cost ($/kWh)','Total cost ($/kWh)'};

    results_firstLaw = [header_firstLaw; num2cell(resultsNumerical_firstLaw)];
    results_secondLaw = [header_secondLaw; num2cell(resultsNumerical_secondLaw)];
    results_dP = [header_dP; num2cell(resultsNumerical_dP)];    
    results_Energy = [header_Energy; num2cell(resultsNumerical_Energy)]; 
    results_Settings = [header_Settings; num2cell(resultsNumerical_Settings)];
    results_Cost = [header_Cost; num2cell(resultsNumerical_Cost)];
    
    %% Writing results to CSV
    clockString = datestr(clock);
    fileString.firstLaw = [clockString, '1: firstLaw.csv'];
    varname.firstLaw = [fileString.firstLaw];
    cell2csv(varname.firstLaw,results_firstLaw);
    
    fileString.secondLaw = [clockString, '2: secondLaw.csv'];
    varname.secondLaw = [fileString.secondLaw];
    cell2csv(varname.secondLaw,results_secondLaw);
    
    fileString.dP = [clockString, '3: dP.csv'];
    varname.dP = [fileString.dP];
    cell2csv(varname.dP,results_dP);
    
    fileString.Energy = [clockString, '4: Energy.csv'];
    varname.Energy = [fileString.Energy];
    cell2csv(varname.Energy,results_Energy);
    
    fileString.Settings = [clockString, '5: Settings.csv'];
    varname.Settings = [fileString.Settings];
    cell2csv(varname.Settings,results_Settings);
    
    fileString.Cost = [clockString, '6: Cost.csv'];
    varname.Cost = [fileString.Cost];
    cell2csv(varname.Cost,results_Cost);

    T_stor_tot = [T_stor_d, T_stor_c];
    fileString.Tstor = [clockString, z, 'T_stor.csv'];
    varname.Tstor = [fileString.Tstor];
    csvwrite(varname.Tstor,T_stor_tot)

    T_SS_tot = [T_SS_d, T_SS_c];
    fileString.TSS = [clockString, z, 'T_SS.csv'];
    varname.TSS = [fileString.TSS];
    csvwrite(varname.TSS,T_SS_tot)

    T_HTF_tot = [T_HTF_d, T_HTF_c];
    fileString.THTF = [clockString, z, 'T_HTF.csv'];
    varname.THTF = [fileString.THTF];
    csvwrite(varname.THTF,T_HTF_tot)
%     
%% Figures
% figure (1)
% hold on
% plot(m_dot_charge_array, resultsNumerical_secondLaw(:,end))
% plot(m_dot_charge_array, resultsNumerical_secondLaw(:,end-2), '-r')
% legend('Without pressure losses', 'With pressure losses')
% xlabel('Mass flow rate (kg/s)'), ylabel('Exergetic efficiency')
% axis([m_dot_charge_array(1) m_dot_charge_array(end) 0 1])

% % figure (2)
% %  Ex_HT(:,end+1) = Ex_HT(:,end);
% % hold on
% % % plot(time,dEx_HT)
% % plot(time_c,Ex_HT(1,:),'r')
% % plot(time_c,Ex_HT(2,:),'b')
% % xlabel('Time (h)');  ylabel('Ex_HT')
% % legend('100 kg/s','500 kg/s')

% figure (3)
% hold on
%  dEx_HT(:,end+1) = dEx_HT(:,end);
% plot(time_c,dEx_HT(1,:),'r')
% plot(time_c,dEx_HT(2,:),'b')
% xlabel('Time (h)');  ylabel('dEx_HT')
% legend('100 kg/s','500 kg/s')

% figure (60)
% hold on
% plot(time_c,results_T_HTF_in(:,1),'b-')
% plot(time_c,results_T_HTF_in(:,2),'r-')
% plot(time_c,results_T_HTF_out(:,1),'b:')
% plot(time_c,results_T_HTF_out(:,2),'r:')
% xlabel('Time (h)');  ylabel('T_HTF_in')
% % legend('100 kg/s','500 kg/s')

% T_HTF_diff = (results_T_HTF_in - results_T_HTF_out)*dt_c;

% figure (61)
% hold on
% plot(time_c,T_HTF_diff(:,1),'b')
% plot(time_c,T_HTF_diff(:,2),'r')
% xlabel('Time (h)');  ylabel('T_HTF_diff')
% legend('100 kg/s','500 kg/s')

% results_dEx_HT(:,1) = log((results_T_HTF_in(:,1)+273.15)/(results_T_HTF_out(:,1)+273.15))*dt_c;
% results_dEx_HT(:,2) = log((results_T_HTF_in(:,2)+273.15)/(results_T_HTF_out(:,2)+273.15))*dt_c;
% %% figure (5)
% % hold on
% % T1 = find(time>1,1);        %For plotting, T1 = 1 hr. mark (can be changed)
% % T2 = find(time>3,1);        %For plotting, T2 = 5 hr. mark (can be changed)
% % T3 = find(time>6,1);        %For plotting, T3 = 10 hr. mark (can be changed)
% % T4 = find(time>9,1);        %For plotting, T3 = 10 hr. mark (can be changed)
% % subplot(2,2,1), plot(tankLength, T_stor(:,T1),'-r',tankLength,T_HTF(:,T1),'-b',tankLength, 390, '--k')
% % title('t = 1 hr')
% % xlabel('Length (m)'), ylabel('Temperature (^oC)')
% % legend('Stor', 'HTF')
% % subplot(2,2,2), plot(tankLength, T_stor(:,T2),'-r',tankLength,T_HTF(:,T2),'-b',tankLength, 390, '--k')
% % title('t = 3 hr')
% % xlabel('Length (m)'), ylabel('Temperature (^oC)')
% % subplot(2,2,3), plot(tankLength, T_stor(:,T3),'-r',tankLength,T_HTF(:,T3),'-b',tankLength, 390, '--k')
% % title('t = 6 hr')
% % xlabel('Length (m)'), ylabel('Temperature (^oC)')
% % subplot(2,2,4), plot(tankLength, T_stor(:,T4),'-r',tankLength,T_HTF(:,T4),'-b',tankLength, 390, '--k')
% % title('t = 9 hr')
% % xlabel('Length (m)'), ylabel('Temperature (^oC)')

% %% figure (6)
% % hold on
% % x_mid = round(x_final/2);   %Index of the middle of tank
% % subplot(2,1,1), plot(time, T_stor(x_mid,:),'--r',time, T_stor(x_final,:),'-r', time, 390, ':k')             
% % xlabel('Time (hours)'), ylabel('Temperature (^oC)')
% % legend('Halfway', 'Exit')
% % title('Storage fluid transient temperature')
% % axis([0 dischargeTime/3600 200 600])
% % subplot(2,1,2), plot(time,T_HTF(x_mid,:),'--b', time, T_HTF(x_final,:),'-b', time, 390, ':k')  
% % xlabel('Time (hours)'), ylabel('Temperature (^oC)')
% % legend('Halfway', 'Exit')
% % title('HTF transient temperature')
% % axis([0 dischargeTime/3600 200 600])

% %% figure (7)
% % plot(time, Q_dot_turbine,'-b')
% % xlabel('Time (hours)'), ylabel('Energy output (MW)')
% % title('Energy output during discharge cycle')
% % axis([0 dischargeTime/3600 0 60])

% %% Figure 8
% figure (8)
% grid on
% plot(time_d,T_stor_d(x_final_d,:),'-r', time_d,T_HTF_d(x_final_d,:),'-b', time_d,T_stor_avg_d,'--r', time_d,m_dot_b_d,':m')
% xlabel('Time (hours)'), ylabel('Temperature (^oC)')
% legend('T_s_t_o_r (exit)', 'T_H_T_F (exit)', 'T_s_t_o_r_,_a_v_g', 'm_b')
% title('Temperatures during discharge cycle')
% axis([0 dischargeTime/3600 0 600])

% %% Figure 9
% % figure (9)
% % plot(tankLength,T_stor(:,t_final),'-r')
% % xlabel('Length (m)'), ylabel('Temperature (^oC)')
% % legend('T_s_t_o_r (final time)')
% % title('Temperatures during discharge cycle')
% % axis([0 L 0 600])
% % 

% %% Figures showing spatial variation of temperature for stor, SS, and HTF

% figure (11)
% hold all
% plot(tankLength_d, T_stor_d(:,round(t_final_d/5)))
% plot(tankLength_d, T_stor_d(:,round(2*t_final_d/5)))
% plot(tankLength_d, T_stor_d(:,round(3*t_final_d/5)))
% plot(tankLength_d, T_stor_d(:,round(4*t_final_d/5)))
% plot(tankLength_d, T_stor_d(:,round(5*t_final_d/5)))
% grid on
% title('Storage fluid temp. at various times during discharge')
% legend('t/t_d = 0.2','t/t_d = 0.4','t/t_d = 0.6','t/t_d = 0.8','t/t_d = 1.0')
% xlabel('Tank length (m)'), ylabel('Temperature (^oC)')
% axis([0 L 200 500])

% figure (10)
% hold all
% plot(time_d, T_stor_d(round(x_final_d/5),:))
% plot(time_d, T_stor_d(round(2*x_final_d/5),:))
% plot(time_d, T_stor_d(round(3*x_final_d/5),:))
% plot(time_d, T_stor_d(round(4*x_final_d/5),:))
% plot(time_d, T_stor_d(round(5*x_final_d/5),:))
% grid on
% title('Storage fluid temp. at various tank locations')
% legend('x/L = 0.2','x/L = 0.4','x/L = 0.6','x/L = 0.8','x/L = 1.0')
% xlabel('Discharge time (hr)'), ylabel('Temperature (^oC)')
% axis([0 dischargeTime/3600 200 500])

% figure (30)
% hold all
% plot(time_d, T_SS_d(round(x_final_d/5),:))
% plot(time_d, T_SS_d(round(2*x_final_d/5),:))
% plot(time_d, T_SS_d(round(3*x_final_d/5),:))
% plot(time_d, T_SS_d(round(4*x_final_d/5),:))
% plot(time_d, T_SS_d(round(5*x_final_d/5),:))
% title('Wall temp. at various tank locations')
% legend('x/L = 0.2','x/L = 0.4','x/L = 0.6','x/L = 0.8','x/L = 1.0')
% xlabel('Discharge time (hr)'), ylabel('Temperature (^oC)')
% axis([0 dischargeTime/3600 200 500])
% 
% figure (31)
% hold all
% plot(time_d, T_HTF_d(round(x_final_d/5),:))
% plot(time_d, T_HTF_d(round(2*x_final_d/5),:))
% plot(time_d, T_HTF_d(round(3*x_final_d/5),:))
% plot(time_d, T_HTF_d(round(4*x_final_d/5),:))
% plot(time_d, T_HTF_d(round(5*x_final_d/5),:))
% title('HTF temp. at various tank locations')
% legend('x/L = 0.2','x/L = 0.4','x/L = 0.6','x/L = 0.8','x/L = 1.0')
% xlabel('Discharge time (hr)'), ylabel('Temperature (^oC)')
% axis([0 dischargeTime/3600 200 500])

%% Figure 32
% figure (32)
% bar(1:3,[E_HTF_frac,E_SS_frac,E_stor_frac])
% figure11_label = {'HTF','Steel','Storage fluid'};
% ylabel('Fraction of total energy stored')
% title('Distribution of total energy stored')
% set(gca,'XTickLabel',figure11_label)
% axis([0 4 0 1])
% grid on

%% References

%[1]  Kolb, G.J.  ?Evaluation of annual performance of 2-tank and thermocline thermal storage systems for trough plants,? Journal of Solar Energy Engineering, Vol. 133, August 2011.
%[2]  Peng, Ding-Yu, Robinson, Donald B., ?A new two-constant equation of state.?  Ind. Eng. Chem. Fundamen., Vol. 15, 1, pp. 59-64, 1976.