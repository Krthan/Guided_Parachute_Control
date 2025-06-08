%% Parameters of the components
%{
This file contains all the parameters that are required for the parachute
simulation - here we split the system into 4 components: 

1 - canopy
2 - suspension lines
3 - PMA (actuators)
4 - Payload

The conversion of the MOI is done in the function call

%}
global atmos_table canopy_radius_uninflated_R0 canopy_radius_inflated_Rp canopy_shape_ratio_epsilon canopy_cop_zp k11 k33 k44 k15 k66...
    system_mass system_com  K system_Ixx system_Iyy system_Izz
%% Density variation wrt Altitude

atmos_table = readtable('atmos.dat.txt', 'ReadVariableNames', false);
atmos_table.Properties.VariableNames = {'Altitude', 'Temperature', 'Pressure', 'Density'};

%% Canopy
%{
We assume the parachute to be hemispheroidal canopy
%}
canopy_mass_m1                          = 22.8; %canopy mass in Kgs
canopy_com_z1                           = -2.76; %center of mass of hemispheroidal canopy in m
canopy_radius_uninflated_R0             = 9.75; %radius of uninflated canopy in m
canopy_radius_inflated_Rp               = 6.50; %radius of inflated canopy in m {Rp = (2/3) * R0}
canopy_shape_ratio_epsilon              = 0.82;
canopy_cop_zp                           = -2.00;  %center of pressure of the canopy from the base {zp = -3/8 * epsilon * Rp}
k11                                     = 0.5; 
k33                                     = 1.0;
k44                                     = 0.24;
k15                                     = 0.75;
k66                                     = 0;
canopy_Iaa                              = 365.93; %MOI along x axis about its own COM
canopy_Ibb                              = 365.93; %MOI along y axis about its own COM
canopy_Icc                              = 623.98; %MOI along z axis about its own COM


%% Suspension lines
suspensionlines_mass_m2                   = 35.3; %suspension lines mass in Kgs
suspensionlines_length_l_SL               = 15.55; %length of suspension lines, l_sl in m
suspensionlines_cone_half_angle_gamma     = 15.31; %cone half angle in degree
suspensionlines_com_z2                    = 7.50; %centroid of suspension lines, z2 in m
suspensionlines_Iaa                       = 1047.09;  %MOI along x axis about its own COM
suspensionlines_Ibb                       = 1047.09;  %MOI along y axis about its own COM
suspensionlines_Icc                       = 770.76;  %MOI along x axis about its own COM



%% PMAs (Pnuematic muscle actuators)
pma_mass_m3                                = 13.4; %four PMAs mass in Kg
pma_length_pressurized_l_PMA               = 5.80; %length of pressurized PMA in m
pma_length_vented_l_star_PMA               = 7.62; %length of vented PMA in m
pma_com_z3                                 = 17.80; %centroid of PMAs system
pma_Iaa                                    = 53.43; %MOI along x axis about its own COM
pma_Ibb                                    = 53.43; %MOI along y axis about its own COM
pma_Icc                                    = 36.98; %MOI along z axis about its own COM



%% Payload
payload_mass_m4                         = 990; %payload mass in Kg
payload_com_z4                          = 21.19; %centroid of payload in m
payload_Iaa                             = 245.59; %MOI along x, y, z axis about its own COM
payload_Ibb                             = 245.59;
payload_Icc                             = 245.59;

%% Calculating Mass and Moment of inertia components with apparent mass

    system_mass = canopy_mass_m1 + suspensionlines_mass_m2 + pma_mass_m3 + payload_mass_m4;

    system_com = (canopy_mass_m1 * canopy_com_z1 + suspensionlines_mass_m2 * suspensionlines_com_z2 + pma_mass_m3 * pma_com_z3 + payload_mass_m4 * payload_com_z4)/system_mass ;

    K = system_mass * system_com;

    canopy_Ixx = canopy_Iaa + canopy_mass_m1 * canopy_com_z1^2;
    canopy_Iyy = canopy_Ibb + canopy_mass_m1 * canopy_com_z1^2;
    canopy_Izz = canopy_Icc;

    suspensionlines_Ixx = suspensionlines_Iaa + suspensionlines_mass_m2 * suspensionlines_com_z2^2;
    suspensionlines_Iyy = suspensionlines_Ibb + suspensionlines_mass_m2 * suspensionlines_com_z2^2;
    suspensionlines_Izz = suspensionlines_Icc;

    pma_Ixx = pma_Iaa + pma_mass_m3 * pma_com_z3^2;
    pma_Iyy = pma_Ibb + pma_mass_m3 * pma_com_z3^2;
    pma_Izz = pma_Icc;

    payload_Ixx = payload_Iaa + payload_mass_m4 * payload_com_z4^2;
    payload_Iyy = payload_Ibb + payload_mass_m4 * payload_com_z4^2;
    payload_Izz = payload_Icc;

    system_Ixx = canopy_Ixx + suspensionlines_Ixx + pma_Ixx + payload_Ixx;
    system_Iyy = canopy_Iyy + suspensionlines_Iyy + pma_Iyy + payload_Iyy;
    system_Izz = canopy_Izz + suspensionlines_Izz + pma_Izz + payload_Izz;


