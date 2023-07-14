function [ alltimes, allwaveforms, allsampnums ] = Get_Spikes_tetrodeICA( datain, data_time, threshold , SampRate, chan)
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
% 0. Do ICA on all chans
data = fastica(datain,'numofic',size(datain,1)); %rows are vars; cols are samp
skw = skewness(data'); % The more skewed the data the larger the spikes probably are
[skwval,sid] = sort(abs(skw));
data = data(sid,:);
% if skewness is > 1, invert spikes
data(skwval>1,:) = -data(skwval>1,:);
% figure(1), for g = 1:size(data,1), plot(data(g,1:3000) + 200*g*ones(size(data(g,1:3000))),'k'), hold on, end
% 1. get threshold crossings -- both under and over
% have to get new threshold since have done ICA
stdout = 3;
threshold = GetThreshold(stdout,data(1,:));
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
padica = 0; % samples on either side

% which other channels to look at?
prox = proximalchans('64D',chan);
allchans = [chan, prox];
nall = length(allchans);
ndur = bef+aft+1;

waveforms = NaN(numwf,(bef+aft+1)*size(data,1));
times = NaN(numwf,1);
artifact = 300;
for g = 1:numwf
     tmpwfshape = data(1,stpt(g):endpt(g));
     if sum(abs(tmpwfshape) > artifact) < 1
         % find minimum
         [minval, minidx] = min(data(1,stpt(g):endpt(g)));
         % grab wf with min in same place
         centerpoint = stpt(g) + minidx;
         newrange = (centerpoint-bef):(centerpoint+aft);
         if (centerpoint - bef) < 0
             newrange = newrange + abs(centerpoint - bef) + 1;
         elseif (centerpoint + aft) > length(data)
             newrange = newrange - ((centerpoint + aft)-length(data));
         end
     
       allwaves = data(:,newrange);
       allwaves = allwaves';
       waveforms(g,:) = allwaves(:);
       times(g) = data_time(centerpoint);
    end % sum(abs(tmpwf) > artifact) < 0
end % for g = 1:numwf

keepwf = ~isnan(waveforms(:,1));
waveforms = waveforms(keepwf,:);
% test if are getting correct waveforms
% newnumwf = size(waveforms,1);
% randrange = ceil(rand(500,1)*newnumwf);
% figure(1)
% plot(waveforms(randrange,:)','k')

end % sub function end

