classdef np_holder < handle
    
    properties
        data_path
        phy_processed_data
        recording_time
        number_of_units
        units
        event_times
    end
    
    methods
        function obj = np_holder(path,varargin)
            
            % set defaults
            default_quality = ["good"];%changing to add more ? 

            % parse varargin
            p = inputParser;
            addParameter(p,'quality_threshold',default_quality,@isstring);
            parse(p,varargin{:});
            quality_threshold = p.Results.quality_threshold;

%             % load and preprocess data
%             f = waitbar(0,'1','Name','Load Data',...
%                 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            % update waitbar
%             waitbar(20/100,f,'loading data...')
            obj.data_path = path;
            
            obj.load_data(path);
            obj.basic_metrics;
            
%             waitbar(60/100,f,'getting cells...')
            obj.get_units(quality_threshold);
            
%             waitbar(100/100,f,'done')

            % delete waitbar
%             delete(f)

        end
        
        function load_data(self,path)
            
            % load numpy files
            phy_path = path + '/continuous/Neuropix-PXI-100.0';
            self.phy_processed_data.spike_times = readNPY(phy_path + '/spike_times.npy');
            self.phy_processed_data.spike_clust = readNPY(phy_path + '/spike_clusters.npy');
            self.phy_processed_data.clust_qual = tdfread(phy_path + '/cluster_group.tsv');
            self.phy_processed_data.clust_info = tdfread(phy_path + '/cluster_info.tsv');
            self.phy_processed_data.spike_templates = readNPY(phy_path + '/templates.npy');

            % convert spike times to time base
            sampFreq = 1/30000;
            self.phy_processed_data.spike_times = double(self.phy_processed_data.spike_times)*sampFreq;
            
        end
        
        function basic_metrics(self)
            % get total recording time (time of last spike) and number of
            % units (clusters with spikes)
            self.recording_time = max(self.phy_processed_data.spike_times);
            self.number_of_units = numel(unique(self.phy_processed_data.spike_clust));
        end
        
        function get_units(self,quality_threshold)
            
            spike_times = self.phy_processed_data.spike_times;
            spike_clust = self.phy_processed_data.spike_clust;
            clust_info = self.phy_processed_data.clust_info;
           

            % get spike times for each cluster 
            clust_ids = unique(spike_clust);
            
            to_del = false(length(clust_ids),1);
            for i = 1:length(clust_ids)
               
                                
                % get cell cluster id, fr, channel, and depth
                self.units(i).cluster_ID = clust_ids(i);
                self.units(i).fr = clust_info.fr(clust_info.id == clust_ids(i));
                self.units(i).channel = clust_info.ch(clust_info.id == clust_ids(i));
                self.units(i).depth = clust_info.depth(clust_info.id == clust_ids(i));
                

                % get cell quality in phy and kilosort
                [~,idx] = min(abs(clust_info.id-double(clust_ids(i))));
                self.units(i).neurontype = convertCharsToStrings(clust_info.neurontype(idx,:)); %adding in neuron type label
                self.units(i).phy_quality = convertCharsToStrings(clust_info.group(idx,:));
                self.units(i).ks_quality = convertCharsToStrings(clust_info.KSLabel(idx,:));

                % strip any spaces in phy quality descriptions
                self.units(i).phy_quality = strtrim(self.units(i).phy_quality);
                
                % get cell spike times
                % unique assures that there are no doubly counted spikes
                % which can happen if cells were merged 
                self.units(i).spike_times = unique(spike_times(spike_clust == clust_ids(i)));
                
                % determine if quality threshold is met
                if ~any(ismember(quality_threshold,self.units(i).phy_quality))
                    to_del(i) = true;
                end
                
            end
            
            % don't include units if quality threshold is not met
            self.units = self.units(to_del == false);
            


        end

        function bin_spikes(self,varargin)
            
            % set defaults
            default_quality = ["good"];
            default_bin_size = 0.010;

            % parse varargin
            validBinSize = @(x) isnumeric(x) && (x > 0);
            p = inputParser;
            addParameter(p,'quality_threshold',default_quality,@isstring);
            addParameter(p,'bin_size',default_bin_size,validBinSize);
            parse(p,varargin{:});
            
            % waitbar
            f = waitbar(0,'1','Name','Bin spikes',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

            % bin data
            quality_threshold = p.Results.quality_threshold;
            bin_size = p.Results.bin_size;

            % set field name to reflect size of bin in milliseconds
            field_name = "binned"; %+string(bin_size*1000)+"ms";
            
            % get bin edges based on size and length of recording (pad with
            % one bin)
            bin_edges = [0:bin_size:self.recording_time+bin_size]';
            
            for cell = 1:length(self.units)
                if any(ismember(quality_threshold,self.units(cell).phy_quality))
                    
                    % update waitbar
                    waitbar(cell/length(self.units),f,'looping through cells...')
                    
                    % count how many spikes fall in each bin
                    values = self.units(cell).spike_times;
                    N = histcounts(values,bin_edges)';
                    
                    % concatenate with times (last times value is unused)
                    self.units(cell).(field_name) = [bin_edges(1:end-1),N];
                end
            end

            % delete waitbar
            delete(f)
        end
        
        function smooth_binned(self,varargin)
            
            % set defaults
            default_quality = ["good"];
            default_sample_time = 0.010;
            default_kernel_size = 0.020;

            % parse varargin
            validNumInput = @(x) isnumeric(x) && (x > 0);
            p = inputParser;
            addParameter(p,'quality_threshold',default_quality,@isstring);
            addParameter(p,'sample_time',default_sample_time,validNumInput);
            addParameter(p,'kernel_size',default_kernel_size,validNumInput);
            parse(p,varargin{:});

            % make waitbar to denote progress
            f = waitbar(0,'1','Name','Smooth firing rates',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            % bin data
            quality_threshold = p.Results.quality_threshold;
            sample_time = p.Results.sample_time;
            kernel_size = p.Results.kernel_size;
            
            % set field name to reflect size of bin in milliseconds
            field_name = "smoothed";%+string(kernel_size*1000)+"ms_kernel_" + string(sample_time*1000)+"ms_sample";
            
            % want times in 1 ms intervals to start
            times = [0:0.001:self.recording_time]';
            
            % make gaussian kernel (from "MATLAB for neuroscientists")
            sigma = kernel_size; %Standard deviation of the kernel = 20 ms
            edges = [-3*sigma:0.001:3*sigma]; %Time ranges form -3*st. dev. to 3*st. dev.
            kernel = normpdf(edges,0,sigma); %Evaluate the Gaussian kernel
            kernel = kernel*0.001; %Multiply by bin width so the probabilities sum to 1

            for cell = 1:length(self.units)
                if any(ismember(quality_threshold,self.units(cell).phy_quality))

                    % update waitbar
                    waitbar(cell/length(self.units),f,'looping through cells...')
                    
                    % first bin at some small interval (here is is 1 microsecond)
                    bin = floor(self.units(cell).spike_times*1000000)/1000000;

                    % get inverse isi estimate of rate
                    isi = bin(2:end)-bin(1:end-1);
                    isi_inv = 1./isi;
                    isi_inv = [0;isi_inv]; % pad with zero

                    % fill every future unfilled millisecond bin with next inv_isi rate value
                    % meaning rate of every bin between two spikes is the inverse of that isi.
                    vq = interp1(bin,isi_inv,times,'next');

                    % set nans at beginning to 0 
                    vq(isnan(vq))=0;

                    %Convolve spike data with the kernel
                    s=conv(vq,kernel); 
                    center = ceil(length(edges)/2); %Find the index of the kernel center
                    s=s(center:end-center+1); %Trim out the relevant portion of the spike density
                    smoothed_rate = [times,s];

                    % down sample to bin size (default = 10 ms)
                    down_samp_ratio = sample_time/0.001; % currently at 1 ms bin
                    smoothed_rate = downsample(smoothed_rate,down_samp_ratio);

                    self.units(cell).(field_name) = smoothed_rate;
                end
                
            end

            % delete waitbar
            delete(f)

        end
        
        function get_waveforms(self,cells,window,num_spikes,type)

%             % parse varargin
%             validNumInput = @(x) isnumeric(x) && (x > 0);
%             p = inputParser;
%             addParameter(p,'quality_threshold',default_quality,@isstring);
%             addParameter(p,'sample_time',default_sample_time,validNumInput);
%             addParameter(p,'kernel_size',default_kernel_size,validNumInput);
%             parse(p,varargin{:});
% 
%             
%             % bin data
%             quality_threshold = p.Results.quality_threshold;
%             sample_time = p.Results.sample_time;
%             kernel_size = p.Results.kernel_size;
%             
%             
%             
            
            
            % waitbar
            f = waitbar(0,'','Name','Extract waveforms',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            path = self.data_path;
            kilosort_dir = convertStringsToChars(path + "/continuous/Neuropix-PXI-100.0");
            % set parameters 
            gwfparams.dataDir = kilosort_dir;      % KiloSort/Phy output folder
            apD = 'continuous.dat';
            gwfparams.fileName = apD;         % .dat file containing the raw 
            gwfparams.dataType = 'int16';     % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = 384;              % Number of channels that were streamed to disk in .dat file
            gwfparams.wfWin = [-40 41];       % Number of samples before and after spiketime to include in waveform
            gwfparams.nWf = num_spikes;                 % Number of spikes to extract                 
            sp = loadKSdir(kilosort_dir);
            
            inc = 0;
            for cellnum = cells

                % update waitbar
                waitbar(inc/length(cells),f,'looping through cells...')
                
                % get cluster ID of cell
                cluster_id = self.units(cellnum).cluster_ID;
                
                all_windows = nan(size(window,1),33*82);
                for idx = 1:size(window,1)
                    t1 = window(idx,1);
                    t2 = window(idx,2);
                    
                    % Vector of cluster spike times (in samples) same length as .spikeClusters
                    gwfparams.spikeClusters = sp.clu(sp.clu==cluster_id);
                    gwfparams.spikeTimes = sp.st(sp.clu==cluster_id); 
                    % get those spikes that are within temporal window of interest
                    gwfparams.spikeClusters = gwfparams.spikeClusters(gwfparams.spikeTimes>= t1 & gwfparams.spikeTimes<= t2);
                    gwfparams.spikeTimes = gwfparams.spikeTimes(gwfparams.spikeTimes>= t1 & gwfparams.spikeTimes<= t2); % Vector of cluster spike times (in samples) same length as .spikeClusters
                    % convert spike times to samples
                    gwfparams.spikeTimes = ceil(gwfparams.spikeTimes*30000);

                    
                    if isempty(gwfparams.spikeTimes)
                        continue
                    end
    
                    
                    % get waveforms
                    wf = getWaveForms(gwfparams);

                    % get the channels to look at waveforms on
                    best_channel = self.units(cellnum).channel;
                    if best_channel<17 % if at bottom of sampled channels of probe
                        nearest_channels = 1:33;
                    elseif best_channel>size( wf.waveForms,3)-16 % if at top of sampled channels
                        nearest_channels = size( wf.waveForms,3)-32:size( wf.waveForms,3);
                    else % otherwise flanking channels
                        nearest_channels = best_channel-16:best_channel+16;
                    end

                    % get average of all spikes for each channel
                    channel_spikes = squeeze(wf.waveForms(:,:,:,:));
                    channel_avgs = squeeze(nanmean(channel_spikes(:,nearest_channels,:),1));

                    % mean subtract
                    channel_avgs = channel_avgs-mean(channel_avgs,2);
                    
                    if type == "concatenated"
                        A = reshape(channel_avgs',1,numel(channel_avgs));
                        all_windows(idx,:) = A;
                    elseif type == "channels"
                        self.units(cellnum).waveforms = channel_avgs;
                    end
                    
                end
                
                self.units(cellnum).waveforms = nanmean(all_windows,1);
                inc = inc+1;
            end
            
            % delete waitbar
            delete(f)

        end
        
        function get_probe_voltage(self,times)
            
            

%             % parse varargin
%             validNumInput = @(x) isnumeric(x) && (x > 0);
%             p = inputParser;
%             addParameter(p,'quality_threshold',default_quality,@isstring);
%             addParameter(p,'sample_time',default_sample_time,validNumInput);
%             addParameter(p,'kernel_size',default_kernel_size,validNumInput);
%             parse(p,varargin{:});

            % make waitbar to denote progress
            f = waitbar(0,'1','Name','Smooth firing rates',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            % bin data
            quality_threshold = p.Results.quality_threshold;
            sample_time = p.Results.sample_time;
            kernel_size = p.Results.kernel_size;
            
            
            
            
            
            % waitbar
            f = waitbar(0,'','Name','Extract waveforms',...
                'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            
            path = data.data_path;
            kilosort_dir = convertStringsToChars(path + "/continuous/Neuropix-PXI-100.0");
            % set parameters 
            gwfparams.dataDir = kilosort_dir;      % KiloSort/Phy output folder
            apD = 'continuous.dat';
            gwfparams.fileName = apD;         % .dat file containing the raw 
            gwfparams.dataType = 'int16';     % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = 384;              % Number of channels that were streamed to disk in .dat file
            gwfparams.wfWin = [-40 41];       % Number of samples before and after spiketime to include in waveform
            gwfparams.nWf = 0;                 % Number of spikes to extract                 
            sp = loadKSdir(kilosort_dir);
            
            
            
            % Vector of cluster spike times (in samples) same length as .spikeClusters
            gwfparams.spikeClusters = sp.clu(sp.clu==cluster_id);
            gwfparams.spikeTimes = sp.st(sp.clu==cluster_id); 
            % get those spikes that are within temporal window of interest
            gwfparams.spikeClusters = gwfparams.spikeClusters(gwfparams.spikeTimes>= t1 & gwfparams.spikeTimes<= t2);
            gwfparams.spikeTimes = gwfparams.spikeTimes(gwfparams.spikeTimes>= t1 & gwfparams.spikeTimes<= t2); % Vector of cluster spike times (in samples) same length as .spikeClusters
            % convert spike times to samples
            gwfparams.spikeTimes = ceil(gwfparams.spikeTimes*30000);

             
            % get waveforms
            wf = getWaveForms(gwfparams);
            
            
            
            
            
            inc = 0;
            for cellnum = cells

                % update waitbar
                waitbar(inc/length(cells),f,'looping through cells...')
                
                % get cluster ID of cell
                cluster_id = self.units(cellnum).cluster_ID;
                
                all_windows = nan(size(window,1),33*82);
                for idx = 1:size(window,1)
                    t1 = window(idx,1);
                    t2 = window(idx,2);
                    

                    
                    if isempty(gwfparams.spikeTimes)
                        continue
                    end
    


                    % get the channels to look at waveforms on
                    best_channel = self.units(cellnum).channel;
                    if best_channel<17 % if at bottom of sampled channels of probe
                        nearest_channels = 1:33;
                    elseif best_channel>size( wf.waveForms,3)-16 % if at top of sampled channels
                        nearest_channels = size( wf.waveForms,3)-32:size( wf.waveForms,3);
                    else % otherwise flanking channels
                        nearest_channels = best_channel-16:best_channel+16;
                    end

                    % get average of all spikes for each channel
                    channel_spikes = squeeze(wf.waveForms(:,:,:,:));
                    channel_avgs = squeeze(nanmean(channel_spikes(:,nearest_channels,:),1));

                    % mean subtract
                    channel_avgs = channel_avgs-mean(channel_avgs,2);
                    
                    if type == "concatenated"
                        A = reshape(channel_avgs',1,numel(channel_avgs));
                        all_windows(idx,:) = A;
                    elseif type == "channels"
                        self.units(cellnum).waveforms = channel_avgs;
                    end
                    
                end
                
                self.units(cellnum).waveforms = nanmean(all_windows,1);
                inc = inc+1;
            end
            
            % delete waitbar
            delete(f)

        end
        
        function get_event_times(self,stim_dur)
            
            %correct to 'spike_times' timebase by subtracting baseline sample number
            %info found in sync_messages.txt file
            ttl_times_samps = readNPY(self.data_path + '/events/Neuropix-PXI-100.0/TTL_1/timestamps.npy');
            text = string(fileread(self.data_path + '/sync_messages.txt'));
            baseline_samps_txt = regexp(text,'0 start time: (\w*)','tokens');
            baseline_samps = str2double(baseline_samps_txt{:});
            ttl_times_samps_corrected = ttl_times_samps - baseline_samps;
            
            % convert to times
            ttl_times=double(ttl_times_samps_corrected)*(1/30000);
            
%             % get pulse interval
%             pulse_int = 1/stim_freq;
%             
%             % get inter event interval
%             iei = diff(ttl_times);
%             
%             % find which interval was greater than pulse interval
%             a=iei>pulse_int;
%             
%             % pad with one so that ones indicate start of stimulus and are
%             % aligned to event times
%             a = [1;a];
%             
%             % select events
%             self.event_times = ttl_times(a==1);
            
            event_times = [];
            t = 0;
            while t < max(ttl_times)
                
                %get stim time, first ttl state change after first gap
                new_inx = find(ttl_times>t,1);
                t = ttl_times(new_inx);
                
                %collate
                event_times = [event_times; t];
                
                %advance t stim_dur + 0.5ms buffer (can modify if necessary)
                t = t + stim_dur + .0005;
                
            end
            self.event_times = event_times;
            
            
        end

        function plot_to_event(self,varargin)
            
            % set defaults
            default_type = "heatmap";
            default_sub = true;
            default_norm = true;
            default_events = self.event_times;


            % parse varargin
            p = inputParser;
            addParameter(p,'type',default_type,@isstring);
            addParameter(p,'baseline_subtract',default_sub);
            addParameter(p,'norm',default_norm);
            addParameter(p,'events',default_events);

            parse(p,varargin{:});
            type = p.Results.type;
            baseline_sub = p.Results.baseline_subtract;
            norm = p.Results.norm;
            events = p.Results.events;

            % make colormap
            ncol = length(self.units);
            colormap = cbrewer('div', 'Spectral', ncol);
            
            tail_pts = 100;
            
            % make array to write data to
            trials = nan(length(events),tail_pts*2+1);
            %idx = nan(length(events),1);
            
            x = ceil(sqrt(length(self.units)));
            
            m_tot = [];
            for cellnum = 1:length(self.units)
                
                %if any(ismember(quality_threshold,self.units(cell).phy_quality))

                for event = 1:length(events)
                    [~,idx] = min(abs(self.units(cellnum).smoothed(:,1)-events(event,1)));
                    if idx-tail_pts >= 0 && idx+tail_pts<= length(self.units(cellnum).smoothed)
                        trials(event,:) = self.units(cellnum).smoothed(idx-tail_pts:idx+tail_pts,2);
                    end
                end
                
                % plot 
                t = -1:0.01:1;
                m = nanmean(trials,1);
                
                if baseline_sub == true
                    m = m- self.units(cellnum).fr;
                end
                if norm == true
                    m = (normalize(m,'range')-0.5)*2;
                end
                
                
                if type == "individual"
                    sd = std(trials,0,1);
                    num = size(trials,1);
                    SEM = sd/sqrt(num);
                    subplot(x,x,cellnum)
                    hold on
                    poly_x = [t,fliplr(t)];
                    poly_y = [m+SEM,fliplr(m-SEM)];
                    h = fill(poly_x,poly_y,colormap(cellnum,:));
                    set(h,'facealpha',.5,'edgealpha',0)
                    plot(t,m,'color',colormap(cellnum,:))
                elseif type == "heatmap"
                    m_tot = [m_tot;m];
                end
            end
            
            if type == "heatmap"
                % plot heatmap
                y = 1:length(self.units);
                hold on
                imagesc(t,y,m_tot); 
                colorbar
                vert = zeros(size(y));
                plot(vert,y,'--','Color',[1 1 1])  
                % set axes
                xlabel('Time (s)')
                xlim([-1,1])
                ylabel('Cell')
                ylim([0.5,cellnum+0.5])
            end
            
        end

    end % end methods
end