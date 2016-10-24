function physio_interactive(varargin)
% physio_interactive - records physio data from MCC USB device in a
% in a text menu-driven way.   
%
% just call the function, and press h to see a list of options.
%
% typically, recording is as easy as starting, then
% 't' to test 10 seconds of data,
% 'd' to set the recording duration in seconds
% 'b' to start to wait for the scanner trigger (the trigger should start in
%       60s or else it will timeout)
% 'p' or 's' to see the timeseries or spectrum.
% 


% check prereqs
if ~strcmp('PCWIN', computer)
    error('Must use 32-bit Windows version of MATLAB.');
end;
if ~exist('analoginput', 'file')
    error('Must have the MATLAB data acquisition toolbox installed.');
end;

cfg.SampleRate = 1000;
cfg.duration = 390;
cfg.channels = [0, 1, 3];
% the only vals recognized are 'resp' 'o2 sat' and 'mr trigger', affects
% plotting and downsampling if chosen
cfg.chan_interp = {'resp', 'o2 sat', 'mr trigger'};
cfg.session_num = 1;
cfg.trig_immed = 0;  %   = 1 means start collecting data when keyboard is pressed
%                        = 0 means use the fMRI trigger 
cfg.downsample = 1;
cfg.downsample_factor = 10;
cfg.spectrum_topfreq = 4;      % when plotting spectrum, what is the highest freq. 


if cfg.trig_immed
    cfg.trigger_type = 'immediate';
else
    cfg.trigger_type = 'software';
    cfg.trigger_cond = 'rising';    % rising edge passing the trigger_cond_val
    cfg.trigger_cond_val = 2;       % (V)olts, this should be the fMRI trigger channel (ADC 3 on USB-1208 device)
    cfg.trigger_chan = 3;           % corresponds to the index of cfg.channels (not the actual device chan num)
end;


while true
    fprintf('Preparing to capture %d ch at %dHz for %d seconds.\n', ...
        length(cfg.channels), cfg.SampleRate, cfg.duration);
    fprintf('Enter (h) for help.\n');
    in_char = input(' -> ', 's');
    
    switch in_char
        case 'b'
            [data_acq, time_acq] = start_capture(cfg);
            cfg.session_num = cfg.session_num + 1;
        
        case 's'
            if exist('time_acq', 'var')
                plot_spectrum(cfg, time_acq, data_acq);
            else
                disp('no data to plot');
            end;
            
            
        case 'p'
            if exist('time_acq', 'var')
                subplot(2, 1, 1);
                plot(time_acq(1:end-1), data_acq(1:end-1, 1));
                hold on
                plot(time_acq(1:end-1), data_acq(1:end-1, 2), 'r');
                subplot(2, 1, 2);
                plot(time_acq(1:end-1), data_acq(1:end-1, 3));
                title(sprintf('run %d', cfg.session_num));
            else
                disp('no data to plot');
            end;
            
        case 'd'
            in_char = input('Enter seconds -> ', 's');
            [dur, conv_ok] = str2num(in_char);
            if conv_ok && isscalar(dur) && dur > 0 && dur < 1000
                cfg.duration = dur;
            end;
        
        case 'r'
            in_char = input('Enter run number -> ', 's');
            [run, conv_ok] = str2num(in_char);
            if conv_ok && isscalar(run)
                cfg.session_num = run;
                fprintf('Run: %d', run);
            end;
            
        case 'h'
            disp('List of commands:');
            disp('  b -> (b)egin capture');
            disp('  e -> (e)val a matlab command');
            disp('  d -> set (d)uration');
            disp('  h -> display (h)elp');
            disp('  q -> (q)uit');
            disp('  p -> (p)lot timecourse');
            disp('  r -> set (r)un number');
            disp('  s -> plot (s)pectrum');
            disp('  t -> (t)est capture for 10 sec');
            % disp('  t -> set (t)rigger mode');
        
        case 'e'
            disp('Enter a blank line or (q) to escape back to data capture');
            while 1
                try
                    curr_cmd = input(' -> ', 's');
                    if strcmp(curr_cmd, 'q') || numel(curr_cmd) == 0
                        break;
                    end;
                    
                    out = evalc(curr_cmd);
                    disp(out);
                catch me
                    disp(me);
                end;
            end;
            
            
        case '1'
            tr_ons = find(diff(data_acq(:,3)) > 2.5);
            %ntr = input('Please choose number of tr''s: ');
            ntr = 16;

            trs_to_plot = tr_ons(1:ntr:end);
            time_plot = mean(diff(time_acq(trs_to_plot)));
            n_tpoints = round(time_plot/mean(diff(time_acq(1:1000))));

            figure;
            for i = 1:length(trs_to_plot)
                ind_range = (1:10:n_tpoints)+trs_to_plot(i);
                if ind_range(end) > length(time_acq)
                    continue;
                end;



                sm_data = conv(data_acq(ind_range,1), 0.1*ones(10,1), 'valid');
                ind_range = ind_range(5:end-5);

                plot(time_acq(ind_range)-time_acq(ind_range(1)), sm_data);
                hold on;
            end;
            hold off;
            
        case 'q'
            break;
            
        case 't'
            % see:
            % http://www.mathworks.com/help/daq/examples/using-analog-input-triggers.html?prodcode=DA&language=en
            % things we need for trigger:
            % mode: immediate, or software
            cfg2 = cfg;
            cfg2.trig_immed = 1;
            cfg2.trigger_type = 'immediate';
            cfg2.duration = 10;
            cfg2.save_file = 0;
            [data_acq, time_acq] = start_capture(cfg2);
            
            % plot data
            subplot(2, 1, 1);
            plot(time_acq(1:end-1), data_acq(1:end-1, 1));
            hold on
            plot(time_acq(1:end-1), data_acq(1:end-1, 2), 'r');
            subplot(2, 1, 2);
            plot(time_acq(1:end-1), data_acq(1:end-1, 3));
            title(sprintf('run %d', cfg.session_num));
            
        otherwise
            fprintf('Do not recognize input "%s".\n\n');
    end;
end; % while

% capture
function [data_acq, time_acq, cfg] = start_capture(cfg)
    ai = analoginput('mcc', 0);
    
    % add channels
    addchannel(ai, cfg.channels);
    
    % set samplerate
    set(ai, 'SampleRate', cfg.SampleRate);
    sr = get(ai, 'SampleRate');
    if sr ~= cfg.SampleRate
        fprintf('setting sample rate to %f instead of %f', ...
            sr, cfg.SampleRate);
        cfg.SampleRate = sr;
    end;
    
    % set duration
    requiredSamples = floor(sr * cfg.duration);
    set(ai, 'SamplesPerTrigger', requiredSamples);

    % grab info
    cfg.dev_info = evalc('disp(ai)');
    
    % log to file
    logfile = sprintf('run_%02d_%s.dat', cfg.session_num, strrep(datestr(now), ' ', '_')); 
    % set logging mode?
    
    % strategy: init output buffer
    data_ind = 1;
    data_acq = zeros(requiredSamples, length(cfg.channels));
    time_acq = zeros(requiredSamples, 1);
    chunk = 1000;
    
  
    % set triggering conditions
    if ~cfg.trig_immed
        set(ai, 'TriggerType', cfg.trigger_type);
        set(ai, 'TriggerRepeat', 0);
        set(ai, 'TriggerCondition', cfg.trigger_cond);
        set(ai, 'TriggerConditionValue', cfg.trigger_cond_val);
        set(ai, 'TriggerChannel', ai.channel(cfg.trigger_chan));
        set(ai, 'TimeOut', 60);
    end;
    
    % start capture
    start(ai);
    while data_ind < requiredSamples
        if requiredSamples < data_ind + chunk
            [data_acq(data_ind:end, :), time_acq(data_ind:end)] = ...
                getdata(ai, requiredSamples-data_ind+1);
            data_ind = data_ind+chunk;
        else
            [data_acq(data_ind:data_ind+chunk-1, :), time_acq(data_ind:data_ind+chunk-1)] = ...
                getdata(ai, chunk);
            data_ind = data_ind+chunk;
        end;
        fprintf('Time elapsed: %.2f and remaining %.2f\n', ...
            data_ind/cfg.SampleRate, (requiredSamples - data_ind)/cfg.SampleRate);
    end; % data loop
        
    % teardown input device
    stop(ai);
    delete(ai);
    
    % should we downsample data?
    % 'mr trigger' channel will be converted to rising edge onsets
    % and all other channels downsampled with matlab 'decimate'
    % timebase will be divided by downsample_factor
    if cfg.downsample
        [cfg, data_acq, time_acq] = data_downsample(cfg, data_acq, time_acq);
    end;
    
    % save data
    % (no spaces or colons in filename)
    if isfield(cfg, 'save_file') && ~cfg.save_file
        return;
    end;
    matfile = strrep(strrep(logfile, '.dat', '.mat'), ':', '-');
    save(matfile, 'cfg', 'data_acq', 'time_acq');
    
return;

% downsampling function
function [cfg, data_acq2, time_acq2] = data_downsample(cfg, data_acq, time_acq)
    % resample timebase
    time_acq2 = time_acq(1:cfg.downsample_factor:end);
    
    % resample trigger
    trig_chan = strcmpi('mr trigger', cfg.chan_interp);
    if isfield(cfg, 'trigger_cond_val')
        trig_thresh = cfg.trigger_cond_val;
    else
        trig_thresh = (max(data_acq(:, trig_chan)) + min(data_acq(:, trig_chan)))/2;
    end;
    
    trig_bin = double(data_acq(:, trig_chan) > trig_thresh);
    trig_ons = round(find(diff(trig_bin) == 1)/cfg.downsample_factor);
    
    data_acq2 = zeros(length(time_acq2), size(data_acq, 2)); 
    
    % resample data
    for i = 1:size(data_acq, 2)
        if size(data_acq2, 1) == ceil(size(data_acq, 1)./cfg.downsample_factor)
            data_acq2(:, i) = decimate(data_acq(:,i), cfg.downsample_factor);
        else
            temp = decimate(data_acq(:,i), cfg.downsample_factor);
            data_acq2(:, i) = temp(:, 1:end-1);
        end;
    end;

    % turn trigger signal into binary
    data_acq2(:, trig_chan) = zeros(size(time_acq2));
    data_acq2(trig_ons, trig_chan) = max(data_acq(:, trig_chan));
    
return;

% plot spectrum
function plot_spectrum(cfg, time_acq, data_acq)
    %if exist(dpss, 'file')
        
    %else
        % use timebase in time_acq to determine SR.
        sr = (length(time_acq)-1)/(time_acq(end) - time_acq(1));
        
        trig_chan = strcmpi('mr trigger', cfg.chan_interp);
        data_acq = data_acq(:, ~trig_chan);
        
        % fft and chop in half
        data_s = abs(fft(data_acq));
        data_s = data_s(1:ceil(size(data_s, 1)/2), :);
        
        % freq bins-> largest is equiv to Nyquist
        fb = linspace(0, sr/2, size(data_s, 1));
        
        % plot only freq range of interest-> up to 4Hz
        last_ind = find(fb > 4, 1, 'first');
        
        plot(fb(2:last_ind), abs(data_s(2:last_ind)));
    %end;
return;


