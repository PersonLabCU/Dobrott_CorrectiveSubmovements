% get events from NPY files
ttl_times_samps = readNPY(pathx + '/events/Neuropix-PXI-100.0/TTL_1/timestamps.npy');

%correct to 'spike_times' timebase by subtracting baseline sample number
%info found in sync_messages.txt file
text = string(fileread(pathx + '/sync_messages.txt'));
baseline_samps_txt = regexp(text,'0 start time: (\w*)','tokens');
baseline_samps = str2double(baseline_samps_txt{:});
ttl_times_samps_corrected = ttl_times_samps - baseline_samps;

% convert to time from samples
sampFreq = 1/30000;
ttl_times = double(ttl_times_samps_corrected)*sampFreq;

            % get pulse interval
            pulse_int = 1/100;
            
            % get inter event interval
            iei = diff(ttl_times);
            
            % find which interval was greater than pulse interval
            a=iei>pulse_int;
            
            % pad with one so that ones indicate start of stimulus and are
            % aligned to event times
            a = [1;a];
            
            % select events
            NP_eventtimes2 = ttl_times(a==1);

data.get_event_times(.05);
NP_eventtimes = data.event_times;

%%

figure; hold on
scatter(ttl_times,ones(length(ttl_times),1),45,'filled')
scatter(NP_eventtimes,ones(length(NP_eventtimes),1)*0.9,45,'filled')
scatter(NP_eventtimes2,ones(length(NP_eventtimes2),1)*1.1,45,'filled')

