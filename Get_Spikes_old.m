function [ times, waveforms ] = Get_Spikes( data, data_time, threshold , SampRate)
% [ times, waveforms ] = Get_Spikes( data, data_time, thresholds , SampRate)
% function that takes a data channel (data)
% the time associated with that channel (data_time),
% a desired threshold value (threshold), 
% and the sampling rate of the data (SampRate)
% and returns the spike times and 2 ms of waveforms
% output args: times, waveforms
% times = Nspikes x 1
% waveforms = Nspikes x 32 (use 32 sampels to represent each wf)
% USAGE
%{

    [ times, waveforms ] = Get_Spikes( data, data_time, threshold , SampRate)

%}

% 1. get threshold crossings
dunder = data <= threshold;
dover = data >= abs(threshold);

if sum(dunder) > sum(dover)
    used = dunder;
    spikesdown = 1;
else
    used = dover;
    spikesdown = 0;
end
dchange = [0,diff(used)];
% grab 12 samples before and 20 samples after first threshold crossing
allsamp = 1:length(data);
if spikesdown == 1;
    idx = allsamp(dchange == -1);
else
    idx = allsamp(dchange == 1);
end

totdur = SampRate*2E-3; % samples for 2 ms of data
bef = floor(1/3*totdur);
aft = floor(2/3*totdur);
idx = idx(idx < length(data)-aft);
idx = idx(idx > bef);
numwf = length(idx);
stpt = idx - bef;
endpt = idx + aft;

waveforms = NaN(numwf,bef+aft+1);
times = NaN(numwf,1);
artifact = 300;
for g = 1:numwf
%     tmpwf = data_time(stpt(g):endpt(g));
    tmpwfshape = data(stpt(g):endpt(g));
    if sum(abs(tmpwfshape) > artifact) < 1
        if spikesdown
            % find minimum
            [minval, minidx] = min(data(stpt(g):endpt(g)));
            % grab wf with min in same place
            
        else
            % find maximum
            [minval, minidx] = max(data(stpt(g):endpt(g)));
        end
        centerpoint = stpt(g) + minidx;
        newrange = (centerpoint-bef):(centerpoint+aft);
        waveforms(g,:) = data(newrange);
        times(g) = data_time(centerpoint);
    end % sum(abs(tmpwf) > artifact) < 0
end % for g = 1:numwf

keepwf = ~isnan(waveforms(:,1));
waveforms = waveforms(keepwf,:);
% test if are getting correct waveforms
newnumwf = size(waveforms,1);
% randrange = ceil(rand(500,1)*newnumwf);
% plot(waveforms(randrange,:)','k')


end % function end

