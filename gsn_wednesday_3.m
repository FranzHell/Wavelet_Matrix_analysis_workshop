
%% generalized eigendecomposition

A = [3 2; 1 3];
B = [1 1; 4 1];
[eigvecs,eigvals] = eig(A,B);

% define vectors
v1 = [.7 .5]';
v2 = eigvecs(:,2);

% matrix-vector multiplication
v1A = A*v1;
v2A = A*v2;
v2B = B*v2;
v2AB = inv(B)*A*v2;

% maximum value for plotting
xval = max([ abs(v1A); abs(v2A) ])*1.1;


figure(1), clf

subplot(221), imagesc(A), axis square, title('A')
subplot(222), imagesc(B), axis square, title('B')


subplot(234)
plot([0 v2(1)],[0 v2(2)],'k','linew',4), hold on
plot([0 v2A(1)],[0 v2A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval]), plot(get(gca,'xlim'),[0 0],'k:'), plot([0 0],get(gca,'ylim'),'k:')
title('Av')


subplot(235)
plot([0 v2(1)],[0 v2(2)],'k','linew',4), hold on
plot([0 v2B(1)],[0 v2B(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval]), plot(get(gca,'xlim'),[0 0],'k:'), plot([0 0],get(gca,'ylim'),'k:')
title('Bv')


subplot(236)
plot([0 v2(1)],[0 v2(2)],'k','linew',4), hold on
plot([0 v2AB(1)],[0 v2AB(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval]), plot(get(gca,'xlim'),[0 0],'k:'), plot([0 0],get(gca,'ylim'),'k:')
title('B^-^1Av')

%% load data

load SSVEP_data

ALLEEG(1).data = double(ALLEEG(1).data);
ALLEEG(2).data = double(ALLEEG(2).data);

%% compute static SSVEP

time2start = 0; % in ms
time2stop  = 1600; % in ms

% specify number of FFT points, which determines the frequency resolution.
% In this case, we can decide the N based on the desired resolution.
% Because the flicker rates are 10 and 12.5 Hz, we need at least .5 Hz resolution.
nFFT = ceil( ALLEEG(1).srate/.25 );

% convert times in ms to indices
[~,timeidx(1)] = min(abs(ALLEEG(1).times-time2start));
[~,timeidx(2)] = min(abs(ALLEEG(1).times-time2stop));
n = diff(timeidx)+1;

% initialize matrices
pow1 = zeros(ALLEEG(1).nbchan,nFFT);
pow2 = zeros(ALLEEG(1).nbchan,nFFT);

% now compute FFT
for chani = 1:ALLEEG(1).nbchan
    fft1 = fft( squeeze(ALLEEG(1).data(chani,timeidx(1):timeidx(2),:)),nFFT )/n;
    fft2 = fft( squeeze(ALLEEG(2).data(chani,timeidx(1):timeidx(2),:)),nFFT )/n;

    % now average FFT power over trials
    pow1(chani,:) = mean( abs( fft1 ).^2 ,2);
    pow2(chani,:) = mean( abs( fft2 ).^2 ,2);
end


% finally, we need the vector of frequencies in Hz
hz = linspace(0,ALLEEG(1).srate/2,floor(nFFT/2)+1);

% plot static SSVEPs from 10 & 12.5 Hz using Headplot function
%frex = [10 12.5];
frex = [6 18];
freqidx = zeros(1,length(frex));
for i=1:2
    [~,freqidx(i)] = min(abs(hz - frex(i)));
end

figure(1), clf
subplot(231)
topoplotIndie(squeeze(pow1(:,freqidx(1))),ALLEEG(1).chanlocs);
set(gca,'clim',[-200 200])
title('10 Hz target')

subplot(232)
topoplotIndie(squeeze(pow2(:,freqidx(1))),ALLEEG(1).chanlocs);
set(gca,'clim',[-80 80])
title('10 Hz flankers')

subplot(233)
topoplotIndie(squeeze((mean(pow1(:,freqidx),2) + mean(pow2(:,freqidx),2))./2),ALLEEG(1).chanlocs);
set(gca,'clim',[-200 200])
title('Average')

%% visually inspect and pick channel for further analysis

chan2use = 'Oz';
chanidx  = strcmpi(chan2use,{ALLEEG(1).chanlocs.labels});

% plot frequency spectrum of EEG for one channel
subplot(212)
plot(hz,pow1(chanidx,1:length(hz)),'.-')
hold on
plot(hz,pow2(chanidx,1:length(hz)),'r.-')
set(gca,'xlim',[1 51])
title([ 'Static SSVEP from channel ' chan2use ])
xlabel('Frequency (Hz)'), ylabel('Power (\muV^2/cm^2)')
legend({'Attention to 10 Hz';'Attention to 12.5 Hz'})

%% nFFT is too short to resolve response at the flicker frequency 

time2start = 0; % in ms
time2stop  = 1600; % in ms 
ntrials = ALLEEG(1).trials; 

% convert times in ms to indices
[~,timeidx(1)] = min(abs(ALLEEG(1).times-time2start));
[~,timeidx(2)] = min(abs(ALLEEG(1).times-time2stop));
n = diff(timeidx)+1;

% define frequency resolution
freq_res = [0.25 1];
colorz = 'rb'; 
figure(2), clf
for i=1:length(freq_res)
    
    % Define length of FFT
    nFFT  = ceil( ALLEEG(2).srate/freq_res(i) );
    
    % now compute FFT for one channel
    fft1 = fft( squeeze(ALLEEG(2).data(chanidx,timeidx(1):timeidx(2),1:ntrials)), nFFT )/n;
    
    % now average FFT power over trials
    pow1 = mean( abs( fft1 ).^2 ,2);
    
    % finally, we need the vector of frequencies in Hz
    hz = linspace(0,ALLEEG(1).srate/2,floor(nFFT/2)+1);
    
    
    % plot frequency spectrum of EEG for one channel
    plot(hz,pow1(1:length(hz)),'.-','Color',colorz(i))
    set(gca,'xlim',[1 20])
    hold on
   
    title('Static SSVEP from channel, variable nFFT')
    xlabel('Frequency (Hz)'), ylabel('Power (\muV^2/cm^2)')
end
legend({['Freq resolution: ' num2str(freq_res(1)) ' Hz'];['Freq resolution: ' num2str(freq_res(2)) ' Hz']},'Location','NorthWest')

%% now for GED

figure(3), clf

nFFT  = ceil( ALLEEG(2).srate/freq_res(1) );

% need the vector of frequencies in Hz
hz = linspace(0,ALLEEG(1).srate/2,floor(nFFT/2)+1);

for condi=1:2
    
    % bandpass filter data around SSVEP frequency
    filtdata = filterFGx(ALLEEG(condi).data,ALLEEG(condi).srate,frex(condi),1);
    contdata = reshape(filtdata,ALLEEG(condi).nbchan,[]);
    contdata = bsxfun(@minus,contdata,mean(contdata,2));
    filtcov  = (contdata*contdata')/(size(contdata,2)-1);
    
    filtdata1 = filterFGx(ALLEEG(condi).data,ALLEEG(condi).srate,frex(condi)-2,1);
    filtdata2 = filterFGx(ALLEEG(condi).data,ALLEEG(condi).srate,frex(condi)+2,1);
    filtdatax = (filtdata1+filtdata2)/2;
    contdata = reshape(filtdatax,ALLEEG(condi).nbchan,[]);
    contdata = bsxfun(@minus,contdata,mean(contdata,2));
    filtcov  = (contdata*contdata')/(size(contdata,2)-1);
    
    %{
    % then broadband covariance
    contdata = reshape(ALLEEG(condi).data,ALLEEG(condi).nbchan,[]);
    contdata = bsxfun(@minus,contdata,mean(contdata,2));
    bbcov    = (contdata*contdata')/(size(contdata,2)-1);
    %}
    
    % GED
    [evecs,evals] = eig(filtcov,bbcov);
    [~,sidx] = sort(real(diag(evals)));
    evecs = real(evecs(:,sidx));
    
    subplot(2,3,condi)
    tmpmap = pinv(evecs');
    topoplotIndie(tmpmap(:,end)./max(abs(tmpmap(:,end))),ALLEEG(1).chanlocs);
    set(gca,'clim',[-1 1])
    geddata = reshape( (contdata'*evecs(:,end))',ALLEEG(condi).pnts,ALLEEG(condi).trials);
    
    % now compute FFT for one channel
    fft1 = fft( squeeze(geddata(timeidx(1):timeidx(2),:)), nFFT )/n;
    
    % now average FFT power over trials
    pow1 = mean( abs( fft1 ).^2 ,2);
    pow1 = pow1./pow1(dsearchn(hz',frex(condi)));
    
    subplot(212), hold on
    plot(hz,pow1(1:length(hz)),'.-')
end



set(gca,'xlim',[1 51])
title([ 'Static SSVEP from channel ' chan2use ])
xlabel('Frequency (Hz)'), ylabel('Power (norm.)')
legend({'Attention to 10 Hz';'Attention to 12.5 Hz'})

%%

%% generalized eigendecomposition to optimize a spatial filter
% find a spatial filter that best differentiates activity between correct
% and error trials. Data were collected in a speeded reaction-time task.
% The first 125 trials are correct responses, the next 125 trials are
% incorrect responses.

load CEdata.mat

% Apply a narrow band-pass filter centered at 7 Hz with a FWHM of 3 Hz.
% Notice how the results depend on the filter frequency! Also try
% commenting out this line to see the broadband activity.
EEG.data = filterFGx(EEG.data,EEG.srate,7,3);


% define time ranges for computing the covariance matrix. You can also try
% different ranges to see how the results are affected.
tidx = dsearchn(EEG.times',[0 700]');


% extract and vectorize data, and compute covariance matrices
datCor = reshape(EEG.data(:,tidx(1):tidx(2),1:125),EEG.nbchan,[]);
datCor = bsxfun(@minus,datCor,mean(datCor,2));
covCor = (datCor*datCor') / size(datCor,2);

datErr = reshape(EEG.data(:,tidx(1):tidx(2),126:250),EEG.nbchan,[]);
datErr = bsxfun(@minus,datErr,mean(datErr,2));
covErr = (datErr*datErr') / size(datErr,2);


% GED
[V,D] = eig(covErr,covCor);

% find eigenvalue order
[~,eigidx] = sort(diag(D));

% compute spatial 'activation pattern'
maps = inv(V');


% Plot the power of the filtered data
figure(4), clf

subplot(221)
topoplotIndie(maps(:,eigidx(end)),EEG.chanlocs);

subplot(222)
% filter is eigenvector with largest eigenvalue. Note that this next line
% also computes the power time series... try deleting abs(hilbert to obtain
% only the real part of the signal.
filter_dataE = abs(hilbert( reshape(EEG.data,EEG.nbchan,[])' * V(:,eigidx(end)) )).^2;
filter_dataE = reshape(filter_dataE,EEG.pnts,EEG.trials);

plot(EEG.times,mean(filter_dataE(:,1:125),2), EEG.times,mean(filter_dataE(:,126:end),2),'linew',2)
xlabel('Time (ms)'), ylabel('Activity (a.u.)')
set(gca,'xlim',[-300 1200])
legend({'Correct response';'Error response'})


subplot(223)
topoplotIndie(maps(:,eigidx(1)),EEG.chanlocs);

subplot(224)
% task+ filter is eigenvector with largest eigenvalue
filter_dataC = abs(hilbert( reshape(EEG.data,EEG.nbchan,[])' * V(:,eigidx(1)) )).^2;
filter_dataC = reshape(filter_dataC,EEG.pnts,EEG.trials);

plot(EEG.times,mean(filter_dataC(:,1:125),2), EEG.times,mean(filter_dataC(:,126:end),2),'linew',2)
xlabel('Time (ms)'), ylabel('Activity (a.u.)')
set(gca,'xlim',[-300 1200])

% finally, adjust y-axis limits
ymin = min([ mean(filter_dataE(:,1:125),2); mean(filter_dataE(:,126:end),2); mean(filter_dataC(:,1:125),2); mean(filter_dataC(:,126:end),2) ]);
ymax = max([ mean(filter_dataE(:,1:125),2); mean(filter_dataE(:,126:end),2); mean(filter_dataC(:,1:125),2); mean(filter_dataC(:,126:end),2) ]);
subplot(222), set(gca,'ylim',[ymin ymax]);
subplot(224), set(gca,'ylim',[ymin ymax]);

%% done.
