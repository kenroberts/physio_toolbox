function varargout = physio_recorder(varargin)
% physio_recorder - records physio data from MCC USB device in a
% scriptable way.  First argument is the command that you want to run. 
%
%
% The easiest way to get started recording data is to just call:
%
% cfg = physio_recorder('getconfig');   % gets a default configuration
% cfg.duration = 180;  % set this to the duration of your run in seconds
% physio_recorder('setconfig', cfg);    % store configuration
% physio_recorder('start');
%
% Commands:
% 
% physio_recorder('getconfig') - returns struct with current config
% information, if 'setconfig' has not been called, this will have
% reasonable defaults
%
% physio_recorder('setconfig', cfg) - sets configuration with values in cfg
% struct
%
% physio_recorder('test') - records 10 secs without waiting for trigger,
% and plots recorded values
%
% physio_recorder('start') - starts recording data in the background
%
% physio_recorder('stop') - saves recording data, and optionally returns
% the data, this function MUST be called.  This function will also block
% until there are enough samples
%
% physio_recorder('trigger') - manually trigger acquisition if it is
% waiting for an fMRI trigger (note: start must have been called)


% if ~strcmp('PCWIN', computer)
%     error('Must use 32-bit Windows version of MATLAB.');
% end;
if ~exist('analoginput', 'file')
    error('Must have the MATLAB data acquisition toolbox installed.');
end;

persistent cfg;

% configure with reasonable defaults
if isempty(cfg)
    cfg = get_config;
end;


cmd_str = varargin{1};

if nargin > 1
    varargin = varargin(2:end);
end;

switch cmd_str
    
    case 'getconfig'
        varargout{1} = cfg;
        
    case 'setconfig'
        if nargin < 2 || ~isstruct(varargin{1})
            error('Must pass in valid config struct as argument');
        end;
        cfg = varargin{1};
        
    case 'test'
        cfg2 = cfg;
        cfg2.trig_immed = 1;
        cfg2.trigger_type = 'immediate';
        cfg2.duration = 10;
        cfg2.save_file = 0;
        cfg2 = start_capture(cfg2);
        [data_acq, time_acq] = stop_capture(cfg2);

        % plot data
        plot_data(cfg2, time_acq, data_acq);
        subplot(2, 1, 1);
        plot(time_acq(1:end-1), data_acq(1:end-1, 1));
        hold on
        plot(time_acq(1:end-1), data_acq(1:end-1, 2), 'r');
        subplot(2, 1, 2);
        plot(time_acq(1:end-1), data_acq(1:end-1, 3));
        title(sprintf('run %d', cfg.session_num));
        
    case 'start'
        cfg = start_capture(cfg);
        
    case 'trigger'
        do_trigger();
    case 'status'
        
    case 'reset'
        
    case 'stop'
        [data_acq, time_acq] = stop_capture(cfg);
        
    otherwise 
        error('Command ''%s'' not recognized.', cmd_str);
end;


return;


% get the current configuration, or reasonable defaults if not set.
function cfg = get_config()

cfg.SampleRate = 1000;
cfg.duration = 390;
cfg.channels = [0, 1, 3];
    % the only vals recognized are 'resp' 'o2 sat' and 'mr trigger', affects
    % plotting and downsampling if chosen
cfg.chan_interp = {'resp', 'o2 sat', 'mr trigger'};
cfg.session_num = 1;
cfg.trig_immed = 0;     %   -> 1 means start collecting data when keyboard is pressed
                        %   -> 0 means use the fMRI trigger 
cfg.downsample = 1;     % 1 -> downsample, 0 -> don't downsample
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

cfg.log_dir = fullfile(getenv('HOMEDRIVE'), getenv('HOMEPATH'), 'MATLAB', 'physio_logs');
if ~exist(cfg.log_dir, 'dir')
    cfg.log_dir = uigetdir('', 'Please choose a directory in which to save log files.');
end;   
 
return;

% start capturing data in background
function cfg = start_capture(cfg)

global ai;
    
ai = analoginput('mcc', 0);

% add channels
addchannel(ai, cfg.channels);

% set samplerate (may update)
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
cfg.logfile = sprintf('run_%02d_%s.mat', cfg.session_num, strrep(datestr(now), ' ', '_'));

% set triggering conditions
if ~cfg.trig_immed
    set(ai, 'TriggerType', cfg.trigger_type);
    set(ai, 'TriggerRepeat', 0);
    set(ai, 'TriggerCondition', cfg.trigger_cond);
    set(ai, 'TriggerConditionValue', cfg.trigger_cond_val);
    set(ai, 'TriggerChannel', ai.channel(cfg.trigger_chan));
    set(ai, 'TimeOut', inf);             % can change timeout
end;

disp('Waiting for trigger...');

% start capture
start(ai);
return;
    
function do_trigger()

global ai;

trigger(ai);

return;


function [data_acq, time_acq] = stop_capture(cfg)
    
global ai

% blocking function
[data_acq, time_acq] = getdata(ai);

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

save(fullfile(cfg.log_dir, cfg.logfile), 'cfg', 'data_acq', 'time_acq');

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


% READ
% https://www.mathworks.com/help/daq/examples/discover-all-other-devices-using-the-legacy-interface.html?prodcode=DA&language=en
% http://www.mathworks.com/help/releases/R2015a/daq/acquire-data-1.html
% https://www.mathworks.com/help/daq/examples/introduction-to-analog-input.html?prodcode=DA&language=en
%
% https://www.mathworks.com/products/daq/code-examples.html

% daqhwinfo('mcc');
% ans = 
%            AdaptorDllName: [1x64 char]
%         AdaptorDllVersion: '2.16 (R2010a)'
%               AdaptorName: 'mcc'
%                BoardNames: {'USB-1208FS'}
%         InstalledBoardIds: {'0'}
%     ObjectConstructorName: {1x3 cell}
    
% ai = analoginput('mcc', 0);
% Display Summary of Analog Input (AI) Object Using 'USB-1208FS'.
%   Acquisition Parameters:  1000 samples per second on each channel.
%                            1000 samples per trigger on each channel.
%                            1 sec. of data to be logged upon START.
%                            Log data to 'Memory' on trigger.
%       Trigger Parameters:  1 'Immediate' trigger(s) on START.
%            Engine status:  Waiting for START.
%                            0 samples acquired since starting.
%                            0 samples available for GETDATA.
% AI object contains no channels.

% http://wiki.biac.duke.edu/biac:experimentalcontrol:biac6hardware
% BioPac Respiratory belt signal connected as MeasComp analog input 0 (ADC0)
% BioPac Cardiac signal connected as MeasComp analog input 1 (ADC1), biopac analog ( o2 ch2, pulse wave ch6, pulse rate ch10 )
% BioPac GSR signal connected as MeasComp analog input 2 (ADC2), biopack analog ch1
% Scanner trigger signal connected as MeasComp analog input 3 (ADC3)

% daqhwinfo(ai) 
% ans = 
%                 AdaptorName: 'mcc'
%                        Bits: 12
%                    Coupling: {'DC Coupled'}
%                  DeviceName: 'USB-1208FS'
%             DifferentialIDs: [0 1 2 3 4 5 6 7]
%                       Gains: []
%                          ID: '0'
%                 InputRanges: [-10 10]
%               MaxSampleRate: 50000
%               MinSampleRate: 0.6000
%              NativeDataType: 'uint16'
%                    Polarity: {'Bipolar'}
%                  SampleType: 'Scanning'
%              SingleEndedIDs: [0 1 2 3 4 5 6 7]
%               SubsystemType: 'AnalogInput'
%               TotalChannels: 8
%     VendorDriverDescription: [1x39 char]
%         VendorDriverVersion: '5'

% props = propinfo(ai);
% props.SampleRate
% ans = 
%                Type: 'double'
%          Constraint: 'bounded'
%     ConstraintValue: [0.6000 50000]
%        DefaultValue: 1000
%            ReadOnly: 'whileRunning'
%      DeviceSpecific: 0

% get HR and respiration
% addchannel(ai, [0 1]);

% get channel desc if needed.
% ch1 = ai.Channel(1);
% get(ch1);

% sampleRate = setverify(ai, 'SampleRate', 1000);

% duration = 5;
% requiredSamples = floor(sampleRate * duration);
% set(ai, 'SamplesPerTrigger', requiredSamples);


% start(ai);
% [d,t] = getdata(ai);
% 
% stop(ai);
% delete(ai);