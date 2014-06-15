function [ m_dot_e, T_powerblock_in, T_powerblock_out, K1, K2, K3, K4, nu ] = powerBlock( powerBlock )

%% Power block specifications [Kolb 2011]

    if strcmp(powerBlock,'config1') == 1   
        m_dot_e = 0.74e6/3600;            %[kg/s] Total mass flow rate to evaporator/power block
        T_powerblock_in = 390;            %[deg C] Design value for power block inlet
        T_powerblock_out = 240;           %[deg C] Design value for power block outlet
        K1 = 0.244444; K2 = 144.66667;    %Coefficients for T_eo curves
        K3 = 0.188889; K4 = -45.66667;    %Coefficients for Q_turbine curves
        nu = 0.37;                        %Power block efficiency
    elseif strcmp(powerBlock,'config2') == 1   
        m_dot_e = 0.98e6/3600;            %[kg/s] Total mass flow rate to evaporator/power block
        T_powerblock_in = 390;            %[deg C] Design value for power block inlet
        T_powerblock_out = 252;           %[deg C] Design value for power block outlet
        K1 = 0.277778; K2 = 143.66667;    %Coefficients for T_eo curves
        K3 = 0.222222; K4 = -53.66667;    %Coefficients for Q_turbine curves
        nu = 0.37;                        %Power block efficiency
    elseif strcmp(powerBlock,'config3') == 1 
        m_dot_e = 1.23e6/3600;            %[kg/s] Total mass flow rate to evaporator/power block
        T_powerblock_in = 390;            %[deg C] Design value for power block inlet
        T_powerblock_out = 263;           %[deg C] Design value for power block outlet
        K1 = 0.311111; K2 = 141.66667;    %Coefficients for T_eo curves
        K3 = 0.244444; K4 = -57.33334;    %Coefficients for Q_turbine curves
        nu = 0.37;                        %Power block efficiency
    elseif strcmp(powerBlock,'config4') == 1 
        m_dot_e = 1.47e6/3600;            %[kg/s] Total mass flow rate to evaporator/power block
        T_powerblock_in = 390;            %[deg C] Design value for power block inlet
        T_powerblock_out = 272;           %[deg C] Design value for power block outlet
        K1 = 0.344444;  K2 = 137.66667;   %Coefficients for T_eo curves
        K3 = 0.288889; K4 = -69.666667;   %Coefficients for Q_turbine curves
        nu = 0.37;                        %Power block efficiency
    elseif strcmp(powerBlock,'config5') == 1 
        m_dot_e = 1.72e6/3600;            %[kg/s] Total mass flow rate to evaporator/power block
        T_powerblock_in = 390;            %[deg C] Design value for power block inlet
        T_powerblock_out = 281;           %[deg C] Design value for power block outlet
        K1 = 0.377778; K2 = 133.66667;    %Coefficients for T_eo curves
        K3 = 0.3; K4 = -71;               %Coefficients for Q_turbine curves
        nu = 0.37;                        %Power block efficiency
    elseif strcmp(powerBlock,'config6') == 1 
        m_dot_e = 1.97e6/3600;              %[kg/s] Total mass flow rate to evaporator/power block
        T_powerblock_in = 390;              %[deg C] Design value for power block inlet
        T_powerblock_out = 289;             %[deg C] Design value for power block outlet
        K1 = 0.433333; K2 = 120.13;         %Coefficients for T_eo curves
        K3 = 0.333334; K4 = -80;            %Coefficients for Q_turbine curves
        nu = 0.37;                          %Power block efficiency
    else
        disp('Select an appropriate power block configuration')
    end                        %Power block efficiency
    
end