function [ filt_data] = butterbandpass( data, n, lowf, highf, samprate )
% implement a nth order butterworth bandpass filter
% [ filt_data] = butterbandpass( data, n, lowf, highf, samprate )
%
% inputs: 
% n = order/2, 
% samprate = amplifier sampling rate
% low f = low cutoff (Hz)
% high f = high cutoff (Hz)
%
% outputs:
% filt_data = filtered data

[A,B,C,D] = butter(n,[lowf highf]/(samprate/2));
sos = ss2sos(A,B,C,D);
g = 1; % gain
filt_data = filtfilt(sos,g,data')';

end

