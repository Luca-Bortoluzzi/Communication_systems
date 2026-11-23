%% Exercise 3 - Part 2: Satellite Link Budget
clc;
clear;

%% 1. System and Geometric Parameters
f0     = 12e9;          % Carrier frequency [Hz] (Ku-band downlink, 12 GHz)
B      = 36e6;          % Transponder bandwidth [Hz] (36 MHz)
PTX    = 50;            % Satellite transmitted power [W]
DTX    = 2;             % Satellite TX antenna diameter [m]
DRX    = 0.6;           % User RX dish diameter [m]
eta_TX = 0.65;          % TX antenna aperture efficiency (typical parabolic reflector)
eta_RX = 0.65;          % RX antenna aperture efficiency (typical residential dish)
Thetae = 40;            % User elevation angle [deg]
Tsys   = 180;           % System noise temperature [K]
Ladd_dB = 3;            % Additional losses (atmospheric, pointing, polarization) [dB]

% Physical constants and orbital geometry
c      = 3e8;           % Speed of light [m/s]
k      = 1.38e-23;      % Boltzmann's constant [J/K]
RE     = 6371e3;        % Earth radius [m]
h      = 35786e3;       % Geostationary (GEO) orbital altitude [m]
lambda = c / f0;        % Signal wavelength [m]

%% 2. Slant Range and Antenna Gains
% Compute satellite-to-user distance d [m] for the selected elevation angle
d = sqrt((RE + h)^2 - RE^2 * cosd(Thetae)^2) - RE * sind(Thetae);

% Compute linear antenna gains: G = eta * (pi * D / lambda)^2
G_TX    = eta_TX * ((pi * DTX) / lambda)^2;
GRX_lin = eta_RX * ((pi * DRX) / lambda)^2;

% Convert antenna gains and TX power to logarithmic units (dBi and dBW)
GTX_dBi = 10 * log10(G_TX);
GRX_dBi = 10 * log10(GRX_lin);
PTX_dBW = 10 * log10(PTX);

%% 3. EIRP and Path Losses
% Effective Isotropic Radiated Power [dBW]
EIRP_dBW = PTX_dBW + GTX_dBi;

% Free-Space Path Loss [dB]: L_FS = 20 * log10(4 * pi * d / lambda)
FSPL_dB = 20 * log10(4 * pi * d / lambda);

% Total link losses including atmospheric, pointing, and polarization margins [dB]
Ltot_dB = Ladd_dB + FSPL_dB;

%% 4. Received Power, Noise Power, and SNR
% Received signal power at the RX antenna output [dBW]
PRX_dBW = EIRP_dBW - Ltot_dB + GRX_dBi;

% Thermal noise power over bandwidth B: P_N = k * T_sys * B [W] -> [dBW]
PN     = k * Tsys * B;
PN_dBW = 10 * log10(PN);

% Signal-to-Noise Ratio at the receiver [dB]
SNR_dB = PRX_dBW - PN_dBW;

%% 5. Link Budget Summary Table
Parameters = {'TX Antenna Gain (GTX)'; 'RX Antenna Gain (GRX)'; ...
    'TX Power (PTX)'; 'EIRP'; 'Free-Space Path Loss (LFS)'; ...
    'Additional Losses (Ladd)'; 'Total Losses (Ltot)'; ...
    'Received Power (PRX)'; 'Noise Power (PN)'; 'SNR'};
Values = [GTX_dBi; GRX_dBi; PTX_dBW; EIRP_dBW; FSPL_dB; ...
    Ladd_dB; Ltot_dB; PRX_dBW; PN_dBW; SNR_dB];
Units = {'dBi'; 'dBi'; 'dBW'; 'dBW'; 'dB'; 'dB'; 'dB'; 'dBW'; 'dBW'; 'dB'};

LinkBudgetTable = table(Parameters, Values, Units);
disp(LinkBudgetTable);