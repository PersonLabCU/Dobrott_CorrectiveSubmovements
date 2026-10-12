% Finds time difference between ttls recorded in NP recording
% (NP_events) and timestamps of ttl generation by motive recording
% system via matlab listeners (ided). 
% Convention: time_dif > 0 means Motive clock started after NP clock
% (typical case)
% The alignment test assumes that NP recorded events include all 
% motive-triggered events (ie if you start motive before the NP recording, 
% this might not work).


function [time_dif] = TimeOffset_NP(NP_events,ided)

% get stim times from motive
Motive_events = ided.m1(ided.m1(:,5) == 1);

% alignment test:
% take diff of each events array.
diffNP = diff(NP_events);
diffM = diff(Motive_events);

%find five longest event intervals in the Motive events array
[longM, longMinx] = maxk(diffM,5);

%loop through Motive event intervals, searching for NP interval with less
%than 1ms error. The reason I'm doing the five top intervals is that if one
%gets messed up in some way (spurious ttl signal for example), then you
%have robustness with the other four intervals. 
time_dif = [];
for n = 1:length(longM)
    test_diff = abs(diffNP - longM(n));
    [min_diff,min_inx] = min(test_diff);
    if min_diff < .001
        time_dif = [time_dif (NP_events(min_inx+1)) - Motive_events(longMinx(n)+1)];
    end
end

%average time_dif for completeness sake. if time_dif is empty, catastrophy
if isempty(time_dif)
    warning('Catastrophic failure. Event arrays do not align.')
else
    time_dif = median(time_dif);
end











