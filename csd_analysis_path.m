
% Code to compute CSD from LFP in response to 2-phase stimuli
% i.e. reversing checkerboard, grating, or full-field
% cmn 06-06

clear; close all; clc

arraytype = '128AN-bottom';

if strcmp(arraytype, '128A')
    load('emap128A.mat');
    num_clus = 2;
elseif strcmp(arraytype, '128AN-bottom')
    load('emap128AN-bottom');
    num_clus = 2;
elseif strcmp(arraytype, '256AN-bottom')
    load('emap256AN-bottom.mat');
    num_clus = 4;
else
    error('Array type is not known.')
end

xp = coordinates(2:end-1, 1);
yp = coordinates(2:end-1, 2);
idx = kmeans(xp, num_clus);

%% grab data files (~ 5 minutes of checkerboard recording)
[amplifier_data,t_amplifier,frequency_parameters,board_dig_in_data,t_dig,path] = get_intan_file_data_path;

%% LFP to CSD
% filter raw data
samprate = frequency_parameters(1).amplifier_sample_rate;
minFreq = 5; maxFreq = 300; % set frequency you want to filter
filtdat = NaN(size(amplifier_data)); % declare empty matrix for filtered data
[b,a] = butter(6,maxFreq/(0.5*samprate)); % low pass filter
[z,p,k] = butter(10,minFreq/(0.5*samprate),'high'); % high pass filter
[sos,g] = zp2sos(z,p,k);
nchaninit = size(amplifier_data,1);
fprintf('Filtering data\n')
[b2,a2] = butter(3,2*[59,61]/samprate,'stop');
[b3,a3] = butter(3,2*[179,181]/samprate,'stop');
d60 = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',59, ...
                 'HalfPowerFrequency2',61, 'DesignMethod','butter','SampleRate',samprate);
d120 = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',119, ...
                 'HalfPowerFrequency2',121, 'DesignMethod','butter','SampleRate',samprate);
d180 = designfilt('bandstopiir','FilterOrder',2, 'HalfPowerFrequency1',179, ...
                 'HalfPowerFrequency2',181, 'DesignMethod','butter','SampleRate',samprate);

textprogressbar('Loading and filtering data: ')                        
for q = 1:nchaninit
    textprogressbar(q*100/nchaninit)
    phi = filtfilt(b ,a,amplifier_data(q,:));     % low pass filter
    lfpfilt = filtfilt(sos,g,phi);     % high pass filter    
    buttLoop = filtfilt(d180,filtfilt(d60,filtfilt(d120,lfpfilt))); % bandpass for 60, 180 Hz
    filtdat(q,:) = buttLoop;
end
textprogressbar('Done :)')

%% stimulus data (2nd channel give on/off ranges)
ch = board_dig_in_data(3,:);
% get actual data frames
stpts = diff(ch);
stpt = find(stpts==1)+1;
stpt = stpt(1:end-2);
endpt = stpt + 1*samprate;
nstim = length(stpt);
sampdur = endpt(1)-stpt(1);

% get stim-triggered average response
data = NaN(nchaninit,sampdur);
textprogressbar('Grabbing stim-triggered response: ')                        
for k = 1:nchaninit
    textprogressbar(k*100/nchaninit)
    alldat = NaN(nstim,sampdur);
    for g = 1:nstim
        alldat(g,:) = filtdat(k,stpt(g):(stpt(g)+sampdur-1));
        xnow = (g-1)*sampdur:(g)*sampdur-1;
    end
    data(k,:) = mean(alldat,1);
end
textprogressbar('Done :)')

%% plot 
time = (1:1:sampdur)*(1/samprate);
col_cycle = {'b', 'r', 'g', 'k', 'm', 'c', [0.5,0.5,0], [0.5,0,0], ...
                [0.3,0,0.5], [0,0.5,0.25], [1,0,1], [1,0.5,0]};
c = length(col_cycle);
figure(1)
for k = 1 : nchaninit
    subplot(1, num_clus, idx(k))
    plot(data(k,:)/2+yp(k), 'Color', col_cycle{mod(k,c)+1}); hold on
    text(500+2*yp(k)+40*xp(k),yp(k),num2str(k), 'Color', col_cycle{mod(k,c)+1})
    title (['Shank: ', num2str(ceil(k/64))])
end
xlabel('Samp #')
ylabel('Location on electrode')

%% Identify bad channels
badchans = [];
prompt = 'Are there any bad channels? 1 = yes, 0 = no\n';
existbad = input(prompt);
if existbad
    prompt = 'Which channels are bad?\n';
    b = input(prompt,'s');
    C = strsplit(b);
    badchans = NaN(length(C),1);
    for k = 1:length(C)
        badchans(k) = str2double(C{k});
    end
end

%% Replace bad channels with average of above and below
for ii = 1 : length(badchans)
    baddepth = yp(badchans(ii));
    badidx = idx(badchans(ii));
    
    yshank = yp(idx == idx(badchans(ii)));
    [sortys, sortysid] = sort(yshank,'descend');

    if baddepth == max(yshank) % one of top channels
       top2 = sortysid(1:2)+64*floor(badchans(ii)/64.1);
       data(badchans(ii), :) = data(top2(badchans(ii)~=top2), :); % set equal to other top
    elseif baddepth == min(yshank) % bottom channel
        prevone = sortysid(end-1)+64*floor(badchans(ii)/64.1); % second bottom channel
        data(badchans(ii), :) = data(prevone,:);
    else
        pos = find(sortysid == badchans(ii)-64*floor(badchans(ii)/64.1));
        top_chan = sortysid(pos-1)+64*floor(badchans(ii)/64.1);
        bottom_chan = sortysid(pos+1)+64*floor(badchans(ii)/64.1);
        data(badchans(ii), :) = (data(top_chan, :) + data(bottom_chan,:))/2;
    end
end 

%% Plot after removing bad channels
figure('Name', 'bad channs replaced')
sgtitle(path);
for k = 1 : nchaninit
    subplot(1, num_clus, idx(k))
    plot(data(k,:)/2+yp(k), 'Color', col_cycle{mod(k,c)+1}); hold on
    text(1000+8*yp(k)+40*xp(k),yp(k),num2str(k), 'Color', col_cycle{mod(k,c)+1})
    title (['Shank: ', num2str(ceil(k/64))])
    ylim([-100, 1100])
    set(gca, 'YTick', 50:100:1050, 'YTickLabel', 1000:-100:0)
end

%% Average all channels at the same depth, do CSD, and plot
uniq_idx = unique(idx);
for ii = 1 : length(uniq_idx)
    idx_shank = uniq_idx(ii);
    inds = find(idx == idx_shank);
    
    x_shank = xp(inds);
    y_shank = yp(inds);
    
    yss = sort(unique(y_shank),'descend');
    avgLFP = NaN(length(yss), size(data,2));
    for k = 1 : length(yss)
        grbch = find(y_shank==yss(k))+64*floor(inds(1)/64.1);
        if length(grbch) > 1
            avgLFP(k,:) = mean(data(grbch,:));
        else
            avgLFP(k,:) = data(grbch,:);
        end
    end

    %%% Plot average LFP signals
    chansp = abs(mean(diff(yss)));
    subv = 0.1; % s - duration of csd plot
    figure
    sgtitle(path);
    subplot(1,6,1)
    subval = round(subv*samprate);
    tt = 1/samprate:1/samprate:subval/samprate;
    for jj = 1 : size(avgLFP,1)
        plot(tt, avgLFP(jj, 1:subval)+yss(jj)); hold on
    end
    title(['avg LFPs, Shank: ', num2str(1+floor(inds(1)/64.1))])
    xlabel('time (s)')
    ylabel('Position on array (um)')
    ylim([-50, 1100])
    set(gca, 'YTick', 50:100:1050, 'YTickLabel', 1000:-100:0)

    %%% do CSD analysis
    maxtime = subv;
    csddur = maxtime*samprate; % s
    interpfactor = 10;
    nchan = size(avgLFP, 1);
    x = linspace(0,subv,interpfactor*csddur);
    for sp = 3 : 7
        subplot(1,6,sp-1)
        sitespace = sp-1;
        CSD = (avgLFP(1:(nchan-2*sitespace),:) + avgLFP((2*sitespace+1):nchan,:)...
                - 2*avgLFP((sitespace+1):(nchan-sitespace),:))/(sitespace^2);
        newc = imresize(CSD,[interpfactor*size(CSD,1) interpfactor*size(CSD,2)],'bilinear');  %%% interpolate at 10x sampling
        y = linspace(max(yss)-chansp*sitespace,min(yss)+chansp*sitespace,interpfactor*size(CSD,1));
        imagesc(x, y, flipud(newc(:,1:interpfactor*csddur)))
        colormap(flipud(jet))
        title(sprintf('%d-site spacing', sp-1))
        ylim([-25, 1100])
        set(gca, 'YTick', linspace(25,1000,11), 'YTickLabel', 0:100:1000)
    end
    colorbar
end

