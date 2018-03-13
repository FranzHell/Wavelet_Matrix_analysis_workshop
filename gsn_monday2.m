%% mikexcohen.com
% phase-based connectivity


%% phase-clustering with different variances

% specify parameters
circ_prop = .1; % proportion of the circle to fill
N = 100; % number of "trials"

% generate phase angle distribution
simdata = rand(1,N) * (2*pi) * circ_prop;


% compute ITPC and preferred phase angle
itpc      = abs(mean(exp(1i*simdata)));
prefAngle = angle(mean(exp(1i*simdata)));


% and plot...
figure(7), clf

% as linear histogram
subplot(3,3,4)
hist(simdata,20)
xlabel('Phase angle'), ylabel('Count')
set(gca,'xlim',[0 2*pi])
title([ 'Observed ITPC: ' num2str(itpc) ])

% and as polar distribution
subplot(1,2,2)
polar([zeros(1,N); simdata],[zeros(1,N); ones(1,N)],'k')
hold on
h = polar([0 prefAngle],[0 itpc],'m');
set(h,'linew',3)
title([ 'Observed ITPC: ' num2str(itpc) ])

%% phase synchronization in empirical data

% load data
load v1_laminar

npnts = size(csd,2);
ntrials = size(csd,3);

%% setup some stuff and junk

% channels for connectivity
chan1idx = 1;
chan2idx = 8;


% create complex Morlet wavelet
cent_freq = 8;
time      = -2:1/srate:2-1/srate;
s         = 8/(2*pi*cent_freq);
wavelet   = exp(2*1i*pi*cent_freq.*time) .* exp(-time.^2./(2*s^2));
half_wavN = (length(time)-1)/2;

% FFT parameters
nWave = length(wavelet)
nData = npnts
nConv = nWave + nData -1

% FFT of wavelet
waveletX = fft(wavelet,nConv);
waveletX = waveletX ./ max(waveletX);

% initialize output time-frequency data
phase_data = zeros(2,length(timevec));
real_data  = zeros(2,length(timevec));


% analytic signal of channel 1
data1X = fft(squeeze(csd(chan1idx,:,1)),nConv);
as = ifft(waveletX.*data1X,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(1,:) = angle(as); % extract phase angles
real_data(1,:)  = real(as); % extract the real part (projection onto real axis)

% analytic signal of channel 1
data1X = fft(squeeze(csd(chan2idx,:,1)),nConv);
as = ifft(waveletX.*data1X,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(2,:) = angle(as);
real_data(2,:)  = real(as);

%% setup figure and define plot handles

% open and name figure
figure(1), set(gcf,'Name','Movie magic minimizes the magic.');

% draw the filtered signals
subplot(221)
filterplotH1 = plot(timevec(1),real_data(1,1),'b');
hold on
filterplotH2 = plot(timevec(1),real_data(2,1),'m');
set(gca,'xlim',[timevec(1) timevec(end)],'ylim',[min(real_data(:)) max(real_data(:))])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(cent_freq) ' Hz' ])

% draw the phase angle time series
subplot(222)
phaseanglesH1 = plot(timevec(1),phase_data(1,1),'b');
hold on
phaseanglesH2 = plot(timevec(1),phase_data(2,1),'m');
set(gca,'xlim',[timevec(1) timevec(end)],'ylim',[-pi pi]*1.1)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angles in polar space
subplot(223)
polar2chanH1 = polar([zeros(1,1) phase_data(1,1)]',repmat([0 1],1,1)','b');
hold on
polar2chanH2 = polar([zeros(1,1) phase_data(2,1)]',repmat([0 1],1,1)','m');
title('Phase angles from two channels')

% draw phase angle differences in polar space
subplot(224)
polarAngleDiffH = polar([zeros(1,1) phase_data(2,1)-phase_data(1,1)]',repmat([0 1],1,1)','k');
title('Phase angle differences from two channels')

%% now update plots at each timestep

for ti=1:5:length(timevec)
    
    % update filtered signals
    set(filterplotH1,'XData',timevec(1:ti),'YData',real_data(1,1:ti))
    set(filterplotH2,'XData',timevec(1:ti),'YData',real_data(2,1:ti))
    
    % update cartesian plot of phase angles
    set(phaseanglesH1,'XData',timevec(1:ti),'YData',phase_data(1,1:ti))
    set(phaseanglesH2,'XData',timevec(1:ti),'YData',phase_data(2,1:ti))
    
    subplot(223)
    cla
    polar([zeros(1,ti) phase_data(1,1:ti)]',repmat([0 1],1,ti)','b');
    hold on
    polar([zeros(1,ti) phase_data(2,1:ti)]',repmat([0 1],1,ti)','m');
    
    subplot(224)
    cla
    polar([zeros(1,ti) phase_data(2,1:ti)-phase_data(1,1:ti)]',repmat([0 1],1,ti)','k');
    
    drawnow
end

%% now quantify phase synchronization between the two channels

% phase angle differences
phase_angle_differences = phase_data(2,:)-phase_data(1,:);

% euler representation of angles
euler_phase_differences = exp(1i*phase_angle_differences);

% mean vector (in complex space)
mean_complex_vector = mean(euler_phase_differences);

% length of mean vector (this is the "M" from Me^ik, and is the measure of phase synchronization)
phase_synchronization = abs(mean_complex_vector);

disp([ 'Synchronization between ' num2str(chan1idx) ' and ' num2str(chan2idx) ' is ' num2str(phase_synchronization) '!' ])

% of course, this could all be done on one line:
phase_synchronization = abs(mean(exp(1i*(phase_data(2,:)-phase_data(1,:)))));

% notice that the order of subtraction is meaningless (see below), which means that this measure of synchronization is non-directional!
phase_synchronization_backwards = abs(mean(exp(1i*(phase_data(1,:)-phase_data(2,:)))));


% now plot mean vector
subplot(224)
hold on
h=polar([0 angle(mean_complex_vector)],[0 phase_synchronization]);
set(h,'linewidth',6,'color','g')

%% phase clustering is phase-invariant

figure(2), clf
subplot(221)
polar(repmat(phase_data(2,:)-phase_data(1,:),1,2)',repmat([0 1],1,length(timevec))','k');
title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(phase_data,1)))))) ])

new_phase_data = phase_data;
for i=2:4
    subplot(2,2,i)
    
    % add random phase offset
    new_phase_data(1,:) = new_phase_data(1,:)+rand*pi;
    
    % plot again
    polar(repmat(new_phase_data(2,:)-new_phase_data(1,:)+pi/2,1,2)',repmat([0 1],1,length(timevec))','k');
    title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(new_phase_data,1)))))) ])
end

%%


%% now in human EEG

% load data and pick channels
load sampleEEGdata

%%

chan1 = 'FCz';
chan2 = 'POz';


min_freq =  2;
max_freq = 50;
num_frex = 30;

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
data1X = fft( reshape(EEG.lap(strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);
data2X = fft( reshape(EEG.lap(strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

% initialize output time-frequency data
itpc = zeros(num_frex,EEG.pnts);

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = nCycs(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX./max(waveletX);
    
    % convolution for chan1
    as1 = ifft(waveletX.*data1X,nConv);
    as1 = as1(half_wave+1:end-half_wave);
    as1 = reshape(as1,EEG.pnts,EEG.trials);
    
    % convolution for chan1
    as2 = ifft(waveletX.*data2X,nConv);
    as2 = as2(half_wave+1:end-half_wave);
    as2 = reshape(as2,EEG.pnts,EEG.trials);
    
    % compute power and average over trials abs(mean(exp(1i*(phase_data(2,:)-phase_data(1,:)))))
    
    itpc(fi,:) = abs(mean(exp(1i * (angle(as1) - angle(as2))),2))
    
end


figure(5), clf
contourf(EEG.times,frex,itpc,40,'linecolor','none')
set(gca,'xlim',[-300 1200],'clim',[0 .3])
colormap hot; colorbar

%% laplacian

% compute laplacian, save as new variable
EEG.lap = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);


time2plot = 250; % ms
chan2plot = 'po3';

figure(4), clf
subplot(221)
topoplotIndie(squeeze(mean(EEG.data(:,dsearchn(EEG.times',time2plot),:),3)),EEG.chanlocs);
title([ 'Voltage (' num2str(time2plot) ')' ])

subplot(222)
topoplotIndie(squeeze(mean(EEG.lap(:,dsearchn(EEG.times',time2plot),:),3)),EEG.chanlocs);
title([ 'Laplacian (' num2str(time2plot) ')' ])

subplot(212)
plot(EEG.times,zscore(squeeze(mean(EEG.data(strcmpi({EEG.chanlocs.labels},chan2plot),:,:),3))), EEG.times,zscore(squeeze(mean(EEG.lap(strcmpi({EEG.chanlocs.labels},chan2plot),:,:),3))),'linew',2)
set(gca,'xlim',[-300 1200])
legend({'Voltage';'Laplacian'})
title([ 'ERP from channel ' chan2plot ])
xlabel('Time (ms)'), ylabel('Data (z-score)')

%% 


%% ISPC over time vs. over trials

% FFT parameters
nWave = length(time);
nData = npnts*ntrials;
nConv = nWave+nData-1;


% initialize output time-frequency data
phase_data = zeros(2,npnts,ntrials);


% FFT of wavelet (need to redo FFT because different n_conv)
waveletX = fft(wavelet,nConv);
waveletX = waveletX ./ max(waveletX);


% analytic signal of channel 1
data1X = fft(reshape(csd(chan1idx,:,:),1,[]),nConv);
as = ifft(waveletX.*data1X,nConv);
as = as(half_wavN+1:end-half_wavN);
as = reshape(as,npnts,ntrials);

% collect real and phase data
phase_data(1,:,:) = angle(as);

% analytic signal of channel 1
data1X = fft(reshape(csd(chan2idx,:,:),1,[]),nConv);
as = ifft(waveletX.*data1X,nConv);
as = as(half_wavN+1:end-half_wavN);
as = reshape(as,npnts,ntrials);

% collect real and phase data
phase_data(2,:,:) = angle(as);


figure(3), clf

% ISPC over trials
subplot(211)
ispc_trials = squeeze( abs(mean(exp(1i*diff(phase_data,[],1)),3)) );
plot(timevec,ispc_trials)
set(gca,'xlim',[-.2 1.2])
xlabel('Time (ms)'), ylabel('ISPC')

% ISPC over time
subplot(212)
ispc_time = squeeze( abs(mean(exp(1i*diff(phase_data,[],1)),2)) );
plot(1:ntrials,ispc_time)
xlabel('Trials'), ylabel('ISPC')

%%


%% cluster vs. phase-lag


% number of time steps
n = 100;

t = [17 86]; % time points where polar distributions are drawn (arb)

figure(4), clf

phasedat = rand(n,1)*pi*(pi/10);
connres  = zeros(n,2);

for i=1:n
    
    % compute phase distribution and eulerize
    phasedat = mod(phasedat + pi/20,2*pi);
    cdd = exp(1i*phasedat);
    
    % phase clustering
    connres(i,1) = abs(mean(cdd));
    
    % phase-lag index
    connres(i,2) = abs(mean( sign(imag(cdd)) ));
    
    
    if i==t(1)
        subplot(223)
        polar([zeros(n,1) phasedat]',[zeros(n,1) ones(n,1)]','k')
    elseif i==t(2)
        subplot(224)
        polar([zeros(n,1) phasedat]',[zeros(n,1) ones(n,1)]','k')
    end
end

subplot(211)
plot(1:n,connres)
xlabel('Time (a.u.)'), ylabel('ISPC or dwPLI')
set(gca,'ylim',[0 1.05],'ytick',0:.2:1)
legend({'ISPC';'PLI'})

hold on
plot([t(1) t(1)],get(gca,'ylim'),'k')
plot([t(2) t(2)],get(gca,'ylim'),'k')

%% ISPC and PLI in empirical data

% (note: try with voltage and laplacian!)

% FFT of data (doesn't change on frequency iteration)
data1X = fft( reshape(EEG.data(strcmpi(chan1,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);
data2X = fft( reshape(EEG.data(strcmpi(chan2,{EEG.chanlocs.labels}),:,:),1,nData) ,nConv);

% initialize output time-frequency data
[ispc,pli] = deal( zeros(num_frex,EEG.pnts) );

% loop over frequencies
for fi=1:num_frex
    
    % create wavelet and get its FFT
    s = nCycs(fi)/(2*pi*frex(fi));
    wavelet  = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX./max(waveletX);
    
    % convolution for chan1
    as1 = ifft(waveletX.*data1X,nConv);
    as1 = as1(half_wave+1:end-half_wave);
    as1 = reshape(as1,EEG.pnts,EEG.trials);
    
    % convolution for chan1
    as2 = ifft(waveletX.*data2X,nConv);
    as2 = as2(half_wave+1:end-half_wave);
    as2 = reshape(as2,EEG.pnts,EEG.trials);
    
    % collect eulerized phase angle differences
    cdd = exp(1i*( angle(as1)-angle(as2) ));
    
    % compute ISPC and PLI (and average over trials!)
    ispc(fi,:) = abs(mean(cdd,2))
    pli(fi,:)  = abs(mean(sign(imag(cdd)),2))
    
end


figure(7), clf
subplot(121)
contourf(EEG.times,frex,ispc,40,'linecolor','none')
set(gca,'xlim',[-300 1200],'clim',[0 .4])
colormap hot; colorbar
title([ 'ISPC between channels ' chan1 ' and ' chan2 ])

subplot(122)
contourf(EEG.times,frex,pli,40,'linecolor','none')
set(gca,'xlim',[-300 1200],'clim',[0 .4])
colormap hot; colorbar
title([ 'PLI between channels ' chan1 ' and ' chan2 ])

%% end.
