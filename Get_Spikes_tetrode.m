function [ alltimes, allwaveforms, allsampnums ] = Get_Spikes_tetrode( data, data_time, threshold , SampRate, chan)
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

    [ times, waveforms ] = Get_Spikes( data, data_time, threshold , SampRate, chan)

%}

% 1. get threshold crossings -- both under and over
dunder = data(1,:) <= threshold; % look at threshold crossing for seed
% dover = data(1,:) >= abs(threshold); % look at threshold crossing for seed

% [timesup, waveformsup, sampnumup] = getsub(dover,data,data_time,SampRate, 0, chan); % spikesdown = 0;
[timesdown, waveformsdown, sampnumdown] = getsub(dunder,data,data_time,SampRate, 1, chan); % spikesdown = 1;
% alltimes = [timesup; timesdown];
alltimes = timesdown;
allwaveforms = waveformsdown;
allsampnums = sampnumdown;
% allwaveforms = [-waveformsup; waveformsdown];
% allsampnums = [sampnumup(:); sampnumdown(:)];
end

function [times, waveforms,sampnum] = getsub(used, data, data_time, SampRate, spikesdown,chan)
dchange = [0,diff(used)];
% grab 12 samples before and 20 samples after first threshold crossing
allsamp = 1:length(data);
if spikesdown == 1;
    idx = allsamp(dchange == -1);
else
    idx = allsamp(dchange == 1);
end

totdur = SampRate*1.5E-3; % samples for 2 ms of data
bef = floor(1/3*totdur);
aft = floor(2/3*totdur);
idx = idx(idx < length(data)-aft);
idx = idx(idx > bef);
numwf = length(idx);
ldata = length(data(1,:));
stpt = idx - bef;
sampnum = stpt;
endpt = idx + aft;

% which other channels to look at?
prox = proximalchans('64D',chan);
allchans = [chan, prox];
nall = length(allchans);

waveforms = NaN(numwf,(bef+aft+1)*nall);
times = NaN(numwf,1);
artifact = 300;
for g = 1:numwf
    for k = 1 % only align waveform 1
        %     tmpwf = data_time(stpt(g):endpt(g));
        tmpwfshape = data(k,stpt(g):endpt(g));
        nl = length(tmpwfshape);
        if sum(abs(tmpwfshape) > artifact) < 1
            if spikesdown
                % find minimum
                [minval, minidx] = min(data(k,stpt(g):endpt(g)));
                % grab wf with min in same place
                
            else
                % find maximum
                [minval, minidx] = max(data(k,stpt(g):endpt(g)));
            end
            centerpoint = stpt(g) + minidx;
            if (centerpoint + aft) <= ldata
                if (centerpoint - bef) > 0
                newrange = (centerpoint-bef):(centerpoint+aft);
                else
                    d = abs(centerpoint - bef);
                    newrange = (centerpoint-bef):(centerpoint+aft);
                    newrange = newrange + d + 1;
                end
            else
                d = ldata - centerpoint;
                newrange = (centerpoint-bef):(centerpoint+aft);
                newrange = newrange - (d+1);
            end
            waveforms(g,(k-1)*nl+1:k*nl) = data(k,newrange);
            times(g) = data_time(centerpoint);
        end % sum(abs(tmpwf) > artifact) < 0
    end % for k = 1 % align first waveform
    % now go through and grab rest of waveforms at same time range
    for k = 2:nall
        % newrange defined for k = 1 
        waveforms(g,(k-1)*nl+1:k*nl) = data(k,newrange); 
    end
end % for g = 1:numwf

keepwf = ~isnan(waveforms(:,1));
waveforms = waveforms(keepwf,:);
% test if are getting correct waveforms
newnumwf = size(waveforms,1);
randrange = ceil(rand(500,1)*newnumwf);
figure(1)
plot(waveforms(randrange,:)','k')

end % sub function end

