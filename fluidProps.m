function [ Tc Pc w M k_stor price_stor A3 A2 A1 A0 ] = fluidProps( storFluid )

%% Fluid properties

    %%%Heat transfer fluid properties [1]
    rho_HTF = 1060;                         %[kg/m^3]  Density of HTF
    cp_HTF = 2446;                          %[J/kgK]  Specific heat of HTF
    k_HTF = 0.098;                          %[W/mK] Thermal conductivity of HTF
    alpha_HTF = k_HTF/(rho_HTF*cp_HTF);     %[m^2/s]  Thermal diffusivity of HTF
    nu_HTF = 9.9e-7;                        %[m^2/s]  Kinematic viscosity of HTF
    price_HTF = 3.96;                       %[$/kg]   Price of HTF

    %%%Storage fluid properties  [2]       
    R = 8.314462175;                        % gas constant [=] J/(mol K)

    if strcmp(storFluid,'Naphthalene') == 1        
        % Naphthalene properties
        Tc = 478.4+273.15;      %[K] Critical temperature
        Pc = 4067.6;            %[kPa]  Critical pressure
        w = 0.309;              % Acentric factor
        M = 128.17;             %[kg/kmol] Molecular weight
        k_stor = 0.13;                                          %[W/mK]  Thermal conductivity of storage fluid
        price_stor = 0.33;                                      %[$/kg]  Price of storage fluid
        A3 = 0; A2 = -1.4e-6; A1 = 4.0443e-3; A0 = -15.48e-3;   %Coefficients for ideal-gas enthalpy function
    elseif strcmp(storFluid,'Sulfur') == 1  
        % Sulfur properties
        Tc = 1040.8+273.15;     %[K] Critical temperature
        Pc = 18208;             %[kPa]  Critical pressure
        w = 0.171;              % Acentric factor
        M = 32.066;             %[kg/kmol] Molecular weight
        k_stor = 0.13;      %[W/mK]  Thermal conductivity of storage fluid
        price_stor = 0.18;  %[$/kg]  Price of storage fluid
        A3 = 3.986e-9/M; A2 = -1.628e-5/M; A1 = 2.218e-2/M; A0 = 27.21/M;  %Coefficients for ideal-gas enthalpy function
    elseif strcmp(storFluid,'pXylene') == 1 
        % %p-Xylene properties
        Tc = 343.35+273.15;     %[K] Critical temperature
        Pc = 3343.73;              %[kPa]  Critical pressure
        w = 0.325;              % Acentric factor
        M = 106.17;           %[kg/kmol] Molecular weight
        k_stor = 0.13;      %[W/mK]  Thermal conductivity of storage fluid
        price_stor = 1.50;  %[$/kg]  Price of storage fluid
        A3 = 8.46599e-10; A2 = -3.68121e-6; A1 = 6.09762e-3; A0 = 3.25151e-1;  %Coefficients for ideal-gas enthalpy function (i.e. curve fit to Cp ideal gas function in kJ/kgK)
    else
        fprintf('Please select an appropriate fluid')
    end 