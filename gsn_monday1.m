

%% create a Morlet wavelet

srate = 1000;         % in hz
time  = -2:1/srate:2; % best practice is to have time=0 at the center of the wavelet
frex  = 6.5;          % frequency of wavelet, in Hz

% create complex sine wave
sine_wave = exp( 1i*2*pi*frex.*time );

% create Gaussian window
s = 20 / (2*pi*frex); % this is the standard deviation of the gaussian
gaus_win  = exp( (-time.^2) ./ (2*s^2) );


% now create Morlet wavelet point-wise multiplication
cmw = sine_wave .* gaus_win;


figure(1), clf

subplot(211)
plot(time,real(cmw))
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'Real part of wavelet at ' num2str(frex) ' Hz' ])

subplot(212)
plot(time,imag(cmw))
xlabel('Time (s)'), ylabel('Amplitude')
title([ 'Imaginary part of wavelet at ' num2str(frex) ' Hz' ])


figure(2), clf
plot3(time,real(cmw),imag(cmw),'linew',3)
axis image
xlabel('Time (s)'), ylabel('Real part'), zlabel('Imaginary part')
rotate3d on

%% Morlet wavelet in the frequency domain

cmwX = fft(cmw)/length(cmw);

hz = linspace(0,srate/2,floor(length(cmw)/2)+1);

figure(3), clf
plot(hz,2*abs(cmwX(1:length(hz))),'o-','linew',3,'markerface','k','markersize',8)
xlabel('Frequency (Hz)'), ylabel('Amplitude')
set(gca,'xlim',[0 frex*4])

%%


%% convolution in a small example

% signal
signal = zeros(1,20);
signal(8:15)=1;

kernel = [1 .8 .6 .4 .2];

% matlab's convolution function
matlab_conv_result = conv(signal,kernel,'same');

figure(4), clf
% plot the signal (impulse or boxcar)
subplot(311)
plot(signal,'o-','linew',2,'markerface','g','markersize',9)
set(gca,'ylim',[-.1 1.1],'xlim',[1 20])
title('Signal')

% plot the kernel
subplot(312)
plot(kernel,'o-','linew',2,'markerface','r','markersize',9)
set(gca,'xlim',[1 20],'ylim',[-.1 1.1])
title('Kernel')

% plot the result of convolution
subplot(313)
plot(matlab_conv_result,'o-','linew',2,'markerface','b','markersize',9)
set(gca,'xlim',[1 20],'ylim',[-.1 3.6])
title('Result of convolution')

%% "manual" time-domain convolution

half_kernel = floor(length(kernel)/2);

% EEG data that we'll use for convolution (must be zero-padded).
dat4conv = [ zeros(1,half_kernel) signal zeros(1,half_kernel) ];

% initialize convolution output
conv_res = zeros(1,length(signal)+length(kernel)-1);

% run convolution
for ti=1:length(dat4conv)-length(kernel)
    tempdata = dat4conv(ti:ti+length(kernel)-1);
    
    % compute dot-product (don't forget to flip the kernel backwards!)
    conv_res(ti+half_kernel) = sum(tempdata*kernel(end:-1:1));
end

% cut off edges
conv_res = conv_res(half_kernel+1:end-half_kernel);

figure(5), clf

plot(conv_res,'o-','linew',2,'markerface','g','markersize',9)
hold on
plot(matlab_conv_result,'o-','linew',2,'markerface','r','markersize',3)
legend({'Time-domain convolution';'Matlab conv function'})

%%


%% convolution as frequency-domain multiplication

%% load in V1 data and pick and electrode, trial, etc.

load v1_laminar.mat

chan2use  = 6;
trial2use = 44;

% extract this bit of data for convenience in the rest of this script
data = squeeze(csd( chan2use,:,trial2use) );

%% create a Morlet wavelet

time = -2:1/srate:2-1/srate; % best practice is to have time=0 at the center of the wavelet
frex = 45; % frequency of wavelet, in Hz

% create complex sine wave
sine_wave = exp(1i*2*pi*frex*time);

% create Gaussian window
s = 7 / (2*pi*frex); % this is the standard deviation of the gaussian
gaus_win  = exp( (-time.^2) ./ (2*s^2) )

% now create Morlet wavelet
cmw  = sine_wave .* gaus_win;

%% define convolution parameters

nData = length(data);
nKern = length(cmw);
nConv = nData + nKern - 1;

%% FFTs

% note that the "N" parameter is the length of convolution, NOT the length
% of the original signals! Super-important!


% FFT of wavelet, and amplitude-normalize in the frequency domain
cmwX = fft(cmw,nConv);
cmwX = cmwX ./ max(cmwX);


% FFT of data
dataX = fft(data,nConv);



% now for convolution...
conv_res = dataX.*cmwX;


% compute hz for plotting
hz = linspace(0,srate/2,floor(nConv/2)+1);

%% some plots...

figure(6), clf

% plot power spectrum of data
subplot(311)
plot(hz,2*abs(dataX(1:length(hz))/length(data)))
set(gca,'xlim',[0 frex*2])

% plot power spectrum of wavelet
subplot(312)
plot(hz,abs(cmwX(1:length(hz))))
set(gca,'xlim',[0 frex*2])

% plot power spectrum of convolution result
subplot(313)
plot(hz,2*abs(conv_res(1:length(hz))/length(data)))
set(gca,'xlim',[0 frex*2])


%% now back to the time domain

% note the differences in lengths of these signals
length(conv_res)
length(data)


% take inverse Fourier transform
conv_res_timedomain = ifft(conv_res);


% cut 1/2 of the length of the wavelet from the beginning and from the end
half_wav = floor( length(cmw)/2 )+1;
conv_res_timedomain = conv_res_timedomain(half_wav-1:end-half_wav);


% and plot
figure(7), clf
plot(timevec,data,'k')
hold on
zoom on
plot(timevec,2*real(conv_res_timedomain),'r','linew',2)
set(gca,'xlim',[-.5 1.5])
legend({'LFP data';'convolution-filtered data'})

%%


%% extract power and phase from many frequencies in wavelet convolution

load sampleEEGdata.mat

% wavelet parameters
min_freq =  2;
max_freq = 40;
num_frex = 37;

channel2use = 'pz';


% set range for variable number of wavelet cycles
range_cycles = [ 4 10 ];

% other wavelet parameters
frex  = logspace(log10(min_freq),log10(max_freq),num_frex);
nCycs = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
time  = -2:1/EEG.srate:2;
half_wave = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = EEG.pnts*EEG.trials;
nConv = nWave+nData-1;


% FFT of data (doesn't change on frequency iteration)
dataX = fft( reshape(EEG.data(strcmpi(channel2use,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

% initialize output time-frequency data
tf = zeros(num_frex,EEG.pnts,2);

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = nCycs(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX./max(waveletX);
    
    % run convolution
    as = ifft(waveletX.*dataX,nConv);
    as = as(half_wave+1:end-half_wave);
as2 = reshape(as,EEG.pnts,EEG.trials);
    
    % compute amplitude (envelope) and average over trials
    tf(fi,:,1) = mean(abs(as2),2);
    
    %compute bandpassfiltered data, envelope and phase for quasi-continous
    %data (original data is epoched data which is used as quasi continous
    %here for demonstration.
    
    tfbandpass(fi,:) =real(as);  
     tfenvelope(fi,:) =abs(as);
     tfphase(fi,:) = angle(as);
    
    % compute ITPC
    tf(fi,:,2) = abs(mean(exp( 1i*(angle(as2))),2));
end


figure 
subplot(2,1,1)
plot(squeeze(tfbandpass(4,1:5000)))
hold on
plot(squeeze(tfenvelope(4,1:5000)))
title('Bandpass and Envelope')
xlabel('Time')
ylabel('Amplitude')
subplot(2,1,2)
plot(squeeze(tfphase(4,1:5000)))
title('Phase')
xlabel('Time')
ylabel('Phase')
% plot results
figure(8), clf

subplot(211)
contourf(EEG.times,frex,squeeze(tf(:,:,1)),40,'linecolor','none')
set(gca,'clim',[0 3],'ydir','normal','xlim',[-300 1000])
title('Power')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

subplot(212)
contourf(EEG.times,frex,squeeze(tf(:,:,2)),40,'linecolor','none')
set(gca,'clim',[0 .6],'ydir','normal','xlim',[-300 1000])
title('ITPC')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

%%
