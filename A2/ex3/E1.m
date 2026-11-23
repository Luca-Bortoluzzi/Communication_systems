
clear; 
clc;

f0 = 12*10^9;   
% 12GHz Falls right in the center of the ITU Region 1 Ku-band 
%               downlink allocation dedicated to Broadcasting Satellite Services 
%               (BSS) and Fixed Satellite Services (FSS)

B = 36*10^6;    
% 36MHz Standard channelization width for Ku-band satellite transponders 
% on major European DTH fleets

PTX = 50;       
% 50W Represents a realistic RF power allocation per carrier/transponder 
% from a satellite Traveling Wave Tube Amplifier (TWTA)

DTX = 2;        
% 2m Typical diameter for a shaped-reflector Ku-band spacecraft downlink 
% antenna designed to generate a continental/regional footprint

DRX = 0.6;      
% 0.6m Standard $60\text{ cm}$ parabolic offset dish for residential 
% DTH television reception in Europe.

Thetae = 40;    
% 40° Realistic mid-latitude look angle for a user in Central/Southern Europe

Tsys = 180;     
% 180K Standard clear-sky receiving system noise temperature for a modern 
% residential Ku-band setup.

RE = 6371; % 6371km
h = 35786; % 35786 km
c  = 3e8;               % Speed of light [m/s]
lambda = c / f0;        % Wavelength [m] 
d = sqrt((RE+h)^2 - RE^2*cosd(Thetae)^2) - RE*sind(Thetae);

theta_vec = 10:0.5:90;
d_vec = sqrt((RE + h)^2 - RE^2 * cosd(theta_vec).^2) - RE * sind(theta_vec);
Lfs_dB = 20 * log10(4 * pi * d_vec / lambda);

figure('Name', 'Slant Range and FSPL vs Elevation Angle');

subplot(2, 1, 1);
plot(theta_vec, d_vec / 1e3, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Elevation Angle \theta_e [deg]');
ylabel('Distance d [km]');
title('Satellite-to-User Slant Range d vs. Elevation Angle');


subplot(2, 1, 2);
plot(theta_vec, Lfs_dB, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Elevation Angle \theta_e [deg]');
ylabel('Free-Space Path Loss L_{FS} [dB]');
title('Free-Space Path Loss vs. Elevation Angle (f_0 = 12 GHz)');


% 1. Monotonic Decrease with Elevation Angle:
%    Both the satellite-to-user slant range d and the free-space path loss L_FS
%    decrease monotonically as the elevation angle theta_e increases from 10 deg
%    to 90 deg. At theta_e = 10 deg (near the horizon), the signal travels the
%    longest path (d = 40581.2 km, L_FS = 206.19 dB). At theta_e = 90 deg
%    (sub-satellite zenith), the user is directly beneath the GEO satellite,
%    so d reduces to the orbital altitude h = 35786 km (L_FS = 205.10 dB).
%
% 2. Small Variation Across the Elevation Range:
%    Because the GEO altitude is much larger than the Earth's radius (h >> R_E),
%    the slant range varies by only ~4795 km (~11.8%) across the entire 10-90 deg
%    interval. Due to the logarithmic scaling of L_FS = 20*log10(4*pi*d/lambda),
%    this translates into a total variation of just 1.09 dB in free-space loss.
%    For our selected elevation angle theta_e = 40 deg, d = 37778.3 km and
%    L_FS = 205.57 dB, sitting roughly midway between the two extremes.
%
% 3. Practical System Implications:
%    Although free-space path loss changes by only ~1.1 dB with elevation angle,
%    low elevation angles (theta_e < 20 deg) must still be avoided in practical
%    DTH systems because the signal traverses a much thicker slice of the
%    troposphere (increasing rain attenuation and gaseous absorption) and the
%    receiving antenna picks up higher thermal noise from the warm ground.