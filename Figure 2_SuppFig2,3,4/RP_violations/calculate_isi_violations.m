function [fpRate,nViol] = calculate_isi_violations(spikeTimes)
% calculate ISI violation based on Hill et al. 2011 (How Single Is My Unit)
% spikeTimes in seconds, sorted
%
% refractoryPeriod = 0.0015;  % 1.5 ms (Allen default)
% minISI           = 0;       % censored period (0 if none)
refractoryPeriod=.8/1000;%.0025/2;
N = length(spikeTimes);
T = spikeTimes(end) - spikeTimes(1);
minISI=1/(30000);
spikeTimes = sort(spikeTimes(:));

isi = diff(spikeTimes);
nViol(1) = sum(isi < minISI);
nViol(2) = sum(isi < refractoryPeriod);

violationWindow = refractoryPeriod;%-minISI;

fpRate = (nViol(2) * T) /(2 * N^2 * violationWindow);


LLobets=1-sqrt(1-(nViol(2) * T) /(N^2 * violationWindow));

LLobets2=.5*(1-sqrt(1-(2*nViol(2) * T) /(N^2 * violationWindow)));

end