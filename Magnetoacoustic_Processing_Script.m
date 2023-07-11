
filename='C:\Loop_fUS_2020-05-27@16-08-08\fUS_block_001.bin'; %Processing file pathway
%parameters of aquisition
lockfreq=5;
Fs=500;
halfswitch=0;
ensemble_length=500;
nchannels=128;
%%
fid = fopen(filename, 'r'); % replace by name of your bin file

tmp = fread(fid, 'double');
fclose(fid);
ensemble_length=500;
nchannels=128;

IQ_tmp = reshape(tmp, [],nchannels*2, ensemble_length);
IQ = IQ_tmp(:,1:nchannels,:)+1i*IQ_tmp(:,(nchannels)+1:nchannels*2,:);


figure
imagesc(squeeze(abs(IQ(:,:,1))))
title('Select ROI')
roiselect=roipoly;

figure
imagesc(squeeze(abs(IQ(:,:,1))))
title('Select BG ROI')
bgroiselect=roipoly;
%%
X=bsxfun(@times,log10(abs(IQ)),roiselect);
XBG=bsxfun(@times,log10(abs(IQ)),bgroiselect);
Xangle=bsxfun(@times,angle(IQ),roiselect);
XBGangle=bsxfun(@times,angle(IQ),bgroiselect);

rtrace=squeeze(sum(sum(X,1),2));
bgtrace=squeeze(sum(sum(XBG,1),2));

figure
subplot(2,1,1)
plot(rtrace)
title('Average Signal in Sample Region')
xlabel('frame number')
subplot(2,1,2)
plot(bgtrace)
title('Average Signal in Background Region')
xlabel('frame number')
%%
T = 1/Fs;             % Sampling period       
L = size(IQ,3);             % Length of signal
t = (0:L-1)*T;        % Time vector

Y = fft(rtrace-mean(rtrace));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

figure
subplot(2,1,1)
plot(f,P1./max(P1)) 
title('Single-Sided Amplitude Spectrum of Sample Region')
xlabel('f (Hz)')
ylabel('|P1(f)|')

Y = fft(bgtrace-mean(bgtrace));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

subplot(2,1,2)
plot(f,P1./max(P1)) 
title('Single-Sided Amplitude Spectrum of Background')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%%

rlockmatsin=zeros(67,128,500);
rlockmatcos=zeros(67,128,500);
for k=1:500
rlockmatsin(:,:,k)=ones(67,128)*sin(2*pi*lockfreq*1/Fs*(k-1));
rlockmatcos(:,:,k)=ones(67,128)*cos(2*pi*lockfreq*1/Fs*(k-1));
end
reallocker=sqrt((mean(abs(IQ).*rlockmatsin,3)).^2+(mean(abs(IQ).*rlockmatcos,3)).^2);
realtheta=atan(mean(abs(IQ).*rlockmatsin,3)./mean(abs(IQ).*rlockmatcos,3));
L10ANS=log10(reallocker);
figure
subplot(2,1,1)
imagesc(log10(reallocker))
title('Magnetomotive')
axis equal
axis tight
subplot(2,1,2)
imagesc(log10(mean(abs(IQ),3)))
title('Averaged B mode')
axis equal
axis tight

