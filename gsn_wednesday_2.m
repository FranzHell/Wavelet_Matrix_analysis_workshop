% introduction to linear algebra and matrix analysis

%% least-squares

n = 10;
% generate some data
b = linspace(1,3,10)' + rand(n,1);

A = [ ones(n,1) (1:n)' ];
x = (A'*A)\A'*b;

%% compute the model-predicted data 

% yHat are the data predicted by the model 
% (design matrix scaled by the coefficients)
yHat = x(1)*A(:,1) + x(2)*A(:,2);


figure(1), clf

% important visual check: plot the original data and model-predicted data
subplot(211)
plot(1:n,b,'o','markerface','b','markersize',20)
hold on
plot(1:n,yHat,'rp-','linew',2,'markersize',30,'markerface','k')
title('Correct')

set(gca,'xlim',[0 n+1],'ylim',[1 4.5])
legend({'observed data';'predicted data'})

%% two common mistakes to avoid

% Swapping the coefficients and design matrix columns.
yHat_error1 = x(2)*A(:,1) + x(1)*A(:,2);

% coefficients scale the y-axis data, not the x-axis data
yHat_error2 = x(1)*b + x(2)*b;

subplot(223)
plot(A(:,2),b,'o','markerface','b','markersize',20)
hold on
plot(A(:,2),yHat_error1,'rp-','linew',2,'markersize',30,'markerface','k')
set(gca,'xlim',[0 n+1],'ylim',[1 max(yHat_error1)+1])
title('Error #1')

subplot(224)
plot(A(:,2),b,'o','markerface','b','markersize',20)
hold on
plot(A(:,2),yHat_error2,'rp-','linew',2,'markersize',30,'markerface','k')
set(gca,'xlim',[0 n+1],'ylim',[1 max(yHat_error1)+1])
title('Error #2')

%% least-squares application in real EEG data


%% load EEG data and extract reaction times in ms

load sampleEEGdata.mat

rts = zeros(size(EEG.epoch));

% loop over trials
for ei=1:EEG.trials
    
    % find the index corresponding to time=0, i.e., trial onset
    [~,zeroloc] = min(abs( cell2mat(EEG.epoch(ei).eventlatency) ));
    
    % reaction time is the event after the trial onset
    rts(ei) = EEG.epoch(ei).eventlatency{zeroloc+1};
end

% create design matrix

A = [ ones(EEG.trials,1) rts' ];
%A(:,2) = A(randperm(length(A)),2);
%% define convolution parameters for time-frequency analysis

freqrange  = [2 20]; % extract only these frequencies (in Hz)
numfrex    = 30;     % number of frequencies between lowest and highest


% set up convolution parameters
wavtime = -2:1/EEG.srate:2;
frex    = linspace(freqrange(1),freqrange(2),numfrex);
nData   = EEG.pnts*EEG.trials;
nKern   = length(wavtime);
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;
nCyc    = logspace(log10(4),log10(12),numfrex);

% create wavelets
cmwX = zeros(numfrex,nConv);
for fi=1:numfrex
    
    % create time-domain wavelet
    s   = nCyc(fi) / (2*pi*frex(fi));
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / (2*s.^2) );
    
    % compute fourier coefficients of wavelet and normalize
    cmwX(fi,:) = fft(cmw,nConv);
    cmwX(fi,:) = cmwX(fi,:) ./ max(cmwX(fi,:));
end


% initialize time-frequency output matrix
tf = zeros(numfrex,EEG.pnts);
tf3d = zeros(numfrex,EEG.pnts,EEG.trials);

% compute Fourier coefficients of EEG data (doesn't change over frequency!)
eegX = fft( reshape(EEG.data(47,:,:),1,[]) ,nConv);

% loop over frequencies
for fi=1:numfrex
    
    % second and third steps of convolution
    as = ifft( cmwX(fi,:).*eegX ,nConv );
    
    % cut wavelet back to size of data
    as = as(halfwav+1:end-halfwav);
    as = reshape(as,EEG.pnts,EEG.trials);
    
    % extract power from all trials
    tf3d(fi,:,:) = abs(as).^2;
    
end % end frequency loop

%% now compute correlations

% reshape the 3D matrix to 2D
tf2d = reshape(tf3d,numfrex*EEG.pnts,EEG.trials)';

% the 2D matrix can be used in a single least squares equation
x = (A'*A)\A'*tf2d;
covmat = reshape(x(2,:),numfrex,EEG.pnts); % demushing

%% show the design and data matrices

figure(2), clf

ax1_h = axes;
set(ax1_h,'Position',[.05 .1 .1 .8])
imagesc(A)
set(ax1_h,'xtick',1:2,'xticklabel',{'Int';'RTs'},'ydir','norm')
ylabel('Trials')


ax2_h = axes;
set(ax2_h,'Position',[.25 .1 .7 .8])
imagesc(tf2d)
set(ax2_h,'ydir','norm','clim',[0 20])
ylabel('Trials')
xlabel('Timefrequency')

colormap gray

%% show the results

figure(3), clf

% show time-frequency map of regressors
contourf(EEG.times,frex,covmat,40,'linecolor','none')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
set(gca,'xlim',[-200 1200],'clim',[-.012 .012])

%%


%% eigendecomposition

A = [3 1; 1 2];
[eigvecs,eigvals] = eig(A);

% define vectors
v1 = [.7 -.5]';
v2 = eigvecs(:,1);

% matrix-vector multiplication
v1A = A*v1;
v2A = A*v2;

% maximum value for plotting
xval = max([ abs(v1A); abs(v2A) ])*1.1;


figure(1), clf

subplot(131)
imagesc(A), axis square
title('Matrix A')


subplot(132)
plot([0 v1(1)],[0 v1(2)],'k','linew',4)
hold on
plot([0 v1A(1)],[0 v1A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval])
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')
legend({'v';'Av'})
title('Not an eigenvector!')



subplot(133)
plot([0 v2(1)],[0 v2(2)],'k','linew',4)
hold on
plot([0 v2A(1)],[0 v2A(2)],'r--','linew',2)
axis square, axis([-xval xval -xval xval])
plot(get(gca,'xlim'),[0 0],'k:')
plot([0 0],get(gca,'ylim'),'k:')
legend({'v';'Av'})
title('Yes an eigenvector!')

%% PCA on simulated data

% data
x = [ 1*randn(1000,1) .4*randn(1000,1) ];

% rotation matrix
th = pi/4;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];

% rotate data
y = x*R1;

% PCA of x (original data)
x = bsxfun(@minus,x,mean(x,1));
covmat = (x'*x) / length(x);
[evecsX,evalsX] = eig(covmat);

% PCA of y (correlated data)
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / length(y);
[evecsY,evalsY] = eig(covmat);


figure(2), clf
% plot original data
subplot(121)
plot(y(:,1),y(:,2),'m.','markersize',5)
set(gca,'xlim',[-5 5],'ylim',[-5 5])
hold on
plot(evalsY(1,1)*[0 evecsY(1,1)],evalsY(1,1)*[0 evecsY(2,1)],'k','linew',4)
plot(evalsY(2,2)*[0 evecsY(1,2)],evalsY(2,2)*[0 evecsY(2,2)],'k','linew',4)
xlabel('x-axis'), ylabel('y-axis')
axis square

% compute component scores
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);

subplot(122)
plot(pc2,pc1,'m.')
set(gca,'xlim',[-5 5],'ylim',[-5 5])
xlabel('PC1 axis'), ylabel('PC2 axis')
axis square

%% vectors vs. values

x = rand(20);
covmat = (x'*x) / length(x);
[evecsX,evalsX] = eig(covmat);

figure(3), clf
subplot(121)
imagesc(evecsX)
axis square

subplot(122)
imagesc(evalsX)
set(gca,'clim',[-.2 .2])
axis square


%% ICA vs. PCA

% generate data

% data
x = [ 1*randn(1000,1) .05*randn(1000,1) ];

% rotation matrix
th = -pi/6;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];
th = -pi/3;
R2 = [ cos(th) -sin(th); sin(th) cos(th) ];

% rotate data
y = [ x*R1 ; x*R2 ];


figure(8), clf
subplot(221)
plot(y(:,1),y(:,2),'o')

datarange = max(y(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('X axis'), ylabel('Y axis')
axis square
title('Data in XY space')


% now PCA
y = bsxfun(@minus,y,mean(y,1));
covmat = (y'*y) / length(y);
[evecsY,evalsY] = eig(covmat);

hold on
plot(evalsY(1,1)*[0 evecsY(1,1)],evalsY(1,1)*[0 evecsY(2,1)],'r','linew',4)
plot(evalsY(2,2)*[0 evecsY(1,2)],evalsY(2,2)*[0 evecsY(2,2)],'r','linew',4)


subplot(222)
pc1 = y*evecsY(:,1);
pc2 = y*evecsY(:,2);

plot(pc2,pc1,'ms')
datarange = max([pc1(:); pc2(:)])*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('PC1 axis'), ylabel('PC2 axis')
axis square
title('Data in PC space')





% now ICA
subplot(223)
plot(y(:,1),y(:,2),'o')
datarange = max(y(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])

ivecs = jader(y');
hold on
plot([0 ivecs(1,1)],[0 ivecs(2,1)],'r','linew',4)
plot([0 ivecs(1,2)],[0 ivecs(2,2)],'r','linew',4)
xlabel('X axis'), ylabel('Y axis')
axis square
title('Data in XY space')



subplot(224)
ic_scores = ivecs*y';
plot(ic_scores(1,:),ic_scores(2,:),'ms')
datarange = max(ic_scores(:))*1.2;
set(gca,'xlim',[-datarange datarange],'ylim',[-datarange datarange])
xlabel('IC1 axis'), ylabel('IC2 axis')
axis square
title('Data in IC space')

%% end.
