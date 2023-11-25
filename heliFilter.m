function [IceSurface, CumDist, d1, d2] = heliFilter(LRF, Lat, Lon, Tc, Md)
% Filter the helicopter movement in EM-Bird laser range finder data
%
% USAGE:
%   [IceSurface, CumDist, d1, d2] = heliFilter2(LRF, Lat, Lon, Tc, Md)
%
%	LRF: Vector containing the elevations measured by the laser range finder
%   Lat: Vector of corresponding measurment points latitude.
%   Lon: Vector of corresponding measurment points longitude.
%   Tc: cut-off period for the filter (in m). Default=100
%   md: max distance between two points for the spline interpolation. Default = 10m
%
%	IceSurface: Vector of filtered elevation representing the estimated ice surface
%   CumDist: Vector of the cummulative distance along the flight track
%   d1: first detrending step
%   d2: second detrending step
%   (detrended altimeter data = altimeter - (d1+d2) )
% Coded by Jean Negrel, Norwegian Polar Instiotute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the filter parameters to default values if not passed as arguments
if nargin<3
    Md = 10;
    if nargin<=2
        Tc = 100;
    end
end

% Remove any data above 24 meters, assuming these correspond to the
% helicopter moving towards calibration elevation.
LRF(LRF>24) = NaN;
% Calculate the cumulative distance along the flight track from Lat and Lon
L = length(Lat);
dist=zeros(L,1);
dist(1) = 0;
for n=1:L-1
    dist(n+1)=latlondistm([Lat(n) Lon(n)],[Lat(n+1) Lon(n+1)]);
end
CumDist = cumsum(dist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First step of the filter
% Filter with butterworth to remove the main trend
idx = ~isnan(LRF);
% interpolate data to an average sampling frequency
Fs = 1/mean(diff(CumDist));
dist = min(CumDist):1/Fs:max(CumDist);
alti = interp1(CumDist(idx), LRF(idx), dist, 'linear');
CumDist = dist;
% Apply the Butterworth filter
Fc = 1/Tc; % Cut-off frequency (1/m)
Fc = Fc/(Fs/2); % Normalised cut-off frequency
[b,a] = butter(2,Fc); % Second order low-pass Butterworth filter
d1 = filtfilt(b,a,alti);
% Detrend the LRF data
step1 = alti - d1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second step of the filter
% Detrend the remaining helicopter motion signal with splines
% Calculate the number of measured points corresponding the the sampling frequency
Md = round(Md/Fs);
% Finding the higher point within the given Md length
M = movmax(step1, Md);
% Constraint the start and end of the selected point to the first and last
% measured point. This limits some edges problems.
M([1;end]) = [step1(1);step1(end)];
% Cleaning the max points selection
idx = (M == step1);
M = M(idx);
% Selecting the correspoinding flight distances
CumDist2 = CumDist(idx);
% Interpolate the spline going through all the selected local max
d2 = interp1(CumDist2, M, CumDist, 'spline');
IceSurface = -(step1-d2);
