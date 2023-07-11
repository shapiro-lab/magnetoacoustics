
%% Magnetoacoustic script modified from UFD script 

clear all, format compact

% Ultrafast acquisition parameters & Image parameters
UF.Depth(1)         = 6;      %initial depth [mm]  
UF.Depth(2)         = 8.2;     % final depth [mm]
UF.NbOfBlocs        = 20;    % number of loop
UF.RcvFreq          = 15.625; % emission frequency
UF.Time_loop        = 2;      % [s]

% UF block parameters
numFrames = 500;
dopAngle = [-6:3:6] * pi/180;
global na;
na = length(dopAngle);
nAccum = 1;
dopFrameRate = 500;
sampling_mode = 2;  % 2 = 100%, 4 = 200%

AntiAliasingFilter  = 20;

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 3;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.connector = 1;

% Specify Trans structure array.
Trans.name = 'L22-14vX';
Trans.units = 'mm';
Trans.frequency = 15.625;
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 25; % mfr data sheet lists 30 Volt limit

% conversion in lambda
UF.Lambda = Resource.Parameters.speedOfSound/UF.RcvFreq*1e-3;   % [mm]
P.startDepth =  floor(UF.Depth(1)/UF.Lambda);   % Acquisition depth in wavelengths
P.endDepth =    ceil(UF.Depth(2)/UF.Lambda);   % This should preferrably be a multiple of 128 samples. Is it true? Charlie

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 1]; % [0.4, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.

% Specify Resources.
maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2)); % this will be used in the receive struct
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 128*ceil(maxAcqLength/128*2)*na*sampling_mode; % this size allows for maximum range %modif
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numFrames;    % Nb of frames stored in RcvBuffer.

% interbuffer 2 used to store IQ data
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  %numFrames % one intermediate buffer needed.
Resource.InterBuffer(1).pagesPerFrame = na*numFrames; %na
Resource.InterBuffer(1).rowsPerFrame = PData(1).Size(1); % this size allows for maximum range
Resource.InterBuffer(1).colsPerFrame = PData(1).Size(2);

Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1);
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);

% Specify Transmit waveform structure.  
% TW(1).type = 'envelope';
% TW(1).envNumCycles = 4;
% TW(1).envFrequency = UF.RcvFreq*ones(1,TW(1).envNumCycles);
% TW(1).envPulseWidth = ones(1,4);%[.15 .7 1 1 .7 .15];
TW(1).type = 'parametric';
TW(1).Parameters = [15.625,.67,6,1]; 

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
               
% - Set event TX attributes.
for n = 1:na   % na transmit events
    TX(n).Steer(1) = dopAngle(n);
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts =  512*ones(1,8); % [330 560 780 1010 1023 1023 1023 1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.  

% For Doppler, use narrow bandwidth coefficients (50% BW) centered at
% 15.625; this is a copy of update function's default coef array for 1
% samples per wave
BPFDop = [ -0.00162 +0.00000 +0.00568 +0.00000 -0.01065 +0.00000 +0.01349 ...
           +0.00000 -0.00858 +0.00000 -0.00955 +0.00000 +0.04312 +0.00000 ...
           -0.08841 +0.00000 +0.13550 +0.00000 -0.17130 +0.00000 +0.18463 ];

switch sampling_mode
    case 4
        sampleMode_str = 'NS200BW';
        Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'TGC', 1, ...
                        'InputFilter', BPFDop,...
                        'sampleMode', sampleMode_str, ...
                        'mode', 0), ...
                        1, 2*na*Resource.RcvBuffer(1).numFrames);
    case 2
        sampleMode_str = 'BS100BW';
        Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'TGC', 1, ...
                        'InputFilter', BPFDop,...
                        'sampleMode', sampleMode_str, ...
                        'demodFrequency', Trans.frequency, ...
                        'mode', 0), ...
                        1, 2*na*Resource.RcvBuffer(1).numFrames);
end
                
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na
        if nAccum>1
            % receive set with replace RF
            Receive(na*2*(i-1)+2*j-1).framenum = i;
            Receive(na*2*(i-1)+2*j-1).acqNum = j;
            % receive set with accum RF
            Receive(na*2*(i-1)+2*j).framenum = i;
            Receive(na*2*(i-1)+2*j).acqNum = j;
            Receive(na*2*(i-1)+2*j).mode = 1;
        else
            Receive(na*(i-1)+j).framenum = i;
            Receive(na*(i-1)+j).acqNum = j;
        end
    end
end

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'replace', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',[],...
                   'regionnum', 1), 1, na);

%_________________________________________________________________________%
% Recon/ReconInfo for reconstruct and store the IQData of each acq        %

% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame.
k = 0;
for frame = 1:numFrames
    Recon(frame) = struct('senscutoff', 0.8, ...
        'pdatanum', 1, ...
        'rcvBufFrame',frame, ...
        'IntBufDest', [1,frame], ...
        'ImgBufDest', [0,0], ...
        'newFrameTimeout', UF.Time_loop*1e3,...  % needed for the software to wait for a new frame in the first rcv slot
        'RINums',(k+1:k+na));
    k = k+na;
    
    % Define ReconInfo structures.
    % We need na ReconInfo structures for na steering angles.
    if  nAccum>1
        for acqNum = 1:na
            ReconInfo(na*(frame-1) + acqNum) = struct('mode', 'replaceIQ', ...   % replace IQ data every time
                'txnum', acqNum, ...
                'rcvnum', 2*na*(frame-1) + 2*acqNum-1, ...       % the index of the receive acquisition
                'pagenum', acqNum,...
                'regionnum', 1);
        end
    else
        for acqNum = 1:na
            ReconInfo(na*(frame-1) + acqNum) = struct('mode', 'replaceIQ', ...   % replace IQ data every time
                'txnum', acqNum, ...
                'rcvnum', na*(frame-1) + acqNum, ...       % the index of the receive acquisition
                'pagenum', acqNum,...
                'regionnum', 1);
        end
    end
end

Recon = assignValues2Fields(Recon,[1:length(Recon)], 'IntBufDest', mat2cell(repmat([1 1], length(Recon),1),ones(length(Recon),1)) );
ReconInfo = assignValues2Fields(ReconInfo, [1:na*numFrames], 'pagenum', [1:na*numFrames]);
%_________________________________________________________________________%
%                                                                         %

% Specify Process structure array.
% process for filtering IQ data
n_external_process = 1;
Process(n_external_process).classname = 'External';                   % process structure for 1st Doppler ensemble
Process(n_external_process).method = 'computeSVD';
Process(n_external_process).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,... % process the most recent frame.
                         'dstbuffer','image',...
                         'dstbufnum', 1,...
                         'dstframenum',1};
                     
EF(1).Function = text2cell('%#EF-1');
n_external_process = n_external_process+1;
                     
% process to save IQ data
Process(n_external_process).classname = 'External';                   % process structure for 1st Doppler ensemble
Process(n_external_process).method = 'saveIQData';
Process(n_external_process).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,... % process the most recent frame.
                         'dstbuffer','none'};
n_external_process = n_external_process+1;
                     
EF(2).Function = text2cell('%#EF-2');

% process to display block time
Process(n_external_process).classname = 'External';                   % process structure for 1st Doppler ensemble
Process(n_external_process).method = 'dispBlockTime';
Process(n_external_process).Parameters = {'srcbuffer','none',... 
                         'dstbuffer','none'};
n_external_process = n_external_process+1;
                     
EF(3).Function = text2cell('%#EF-3');

% Specify SeqControl structure arrays. 
maxReturnTrip = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2)*2;
maxTOF = round(2*maxReturnTrip/Trans.frequency);  % to be secure, 2 times
time_compound = nAccum*na*maxTOF;
time_ensemble = 1/(dopFrameRate*1e-6)-(time_compound-maxTOF);
time_2_nextSEQ = UF.Time_loop*1e6- 1/(dopFrameRate*1e-6)*numFrames - time_ensemble;

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start, not used here
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = maxTOF;  % 200 usecs
SeqControl(3).command = 'triggerOut';
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'timeToNextAcq';  % time between frames
SeqControl(5).argument = time_ensemble;  % 2 msec (<=> 500 Hz)
SeqControl(6).command = 'timeToNextAcq';  % time between blocks
SeqControl(6).argument = time_2_nextSEQ;  % 3 sec 
SeqControl(7).command = 'sync';  % synchronisation soft hard, (no ping pong)
SeqControl(7).argument = UF.Time_loop*1e6;  % time out for sync  % needed for the software to wait for the hardware to end the last TTNA of the loop
SeqControl(8).command = 'loopCnt'; % - Set loop count. for looping on blocs
SeqControl(8).argument = UF.NbOfBlocs-1;  %
SeqControl(8).condition = 'counter1';
SeqControl(9).command = 'loopTst';  % loop test
SeqControl(9).argument = [];    % set apres
SeqControl(9).condition = 'counter1';

nsc = 10; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

% set loop count
Event(n).info = 'start counter';
Event(n).tx = 0;   % use next TX structure.
Event(n).rcv = 0;
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = 8;
n = n+1;
SeqControl(9).argument = n;


for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:na                    % Acquire frame
        if nAccum>1
            for k = 1:nAccum
                if k == 1
                    Event(n).info = 'Acquire RF';
                    Event(n).tx = j;   % use next TX structure.
                    Event(n).rcv = 2*na*(i-1)+(j-1)*2+1;
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = [2,3];
                    n = n+1;
                else
                    Event(n).info = 'Acquire RF';
                    Event(n).tx = j;   % use next TX structure.
                    Event(n).rcv = 2*na*(i-1)+(j-1)*2+2;
                    Event(n).recon = 0;      % no reconstruction.
                    Event(n).process = 0;    % no processing
                    Event(n).seqControl = [2];
                    n = n+1;
                end
            end
        else
            Event(n).info = 'Acquire RF';
            Event(n).tx = j;   % use next TX structure.
            Event(n).rcv = na*(i-1)+j;
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            if j==1
                Event(n).seqControl = [2,3]; %Send Trigger at start of each block
            else
                Event(n).seqControl=[2];
            end
            n = n+1;
        end
    end

    % set sync hardware and software for the first TX
    Event(2).seqControl = [7,2,3];

    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [5,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer all acqs to host buffer
      nsc = nsc+1;

    Event(n).info = 'Reconstruct'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = i; % 1;      % reconstruction
    Event(n).process = 0;    % process
    Event(n).seqControl = 0;%4;
    n = n+1;
end

Event(n-2).seqControl = [6,nsc-1]; % modify last acquisition of the lat frame TTNA

Event(n).info = 'Sequence Time Control'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 3;    % process
Event(n).seqControl = 0;
n = n+1;

% Event(n).info = 'External SVD Process'; 
% Event(n).tx = 0;         % no transmit
% Event(n).rcv = 0;        % no rcv
% Event(n).recon = 0;      % no reconstruction
% Event(n).process = 1;    % process
% Event(n).seqControl = 0;
% n = n+1;

Event(n).info = 'save IQ Data Process'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 2;    % process
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'loop everything'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 9 ;
n = n+1;

% Rcv profile
RcvProfile.DCsubtract = 'on';               % substract DC signal if 'on'
RcvProfile.AntiAliasCutoff = AntiAliasingFilter;    % antialiasing analogical filter cuttoff freq [MHz]
RcvProfile.PgaGain = 30;     %24 ou 30              % analog gain in dBfor preamp  %#test
RcvProfile.LnaGain = 24;                   % gain in dB of the fixed gain low noise amp
RcvProfile.LnaZinSel = 31;    % Force high-Z state for best Doppler sensitivity  %#test

%Check coherence of the VS structures
if max([Event(:).rcv])>length(Receive)
    warning('Probably a problem of receive indexing in the Event struct')
elseif max([Event(:).tx])>length(TX)
    warning('Probably a problem of TX indexing in the Event struct')
elseif max([Event(:).seqControl])>length(SeqControl)
    warning('Probably a problem of SeqControl indexing in the Event struct')
elseif max([ReconInfo(:).txnum])>length(TX)
    warning('Probably a problem of TX indexing in the ReconInfo struct')
% elseif max([ReconInfo(:).pagenum])>numFrames
%     warning('Probably a problem of interbuffer pages management in the ReconInfo struct')
elseif max([ReconInfo(:).rcvnum])>length(Receive)
    warning('Probably a problem of receive indexing in the ReconInfo struct')
end

% UI controls

% Set TPCHighVoltage for profile one to 20V
UI(1).Statement = '[result,hv] = setTpcProfileHighVoltage(3,1);';
UI(2).Statement = 'hv1Sldr = findobj(''Tag'',''hv1Sldr'');';
UI(3).Statement = 'set(hv1Sldr,''Value'',hv);';
UI(4).Statement = 'hv1Value = findobj(''Tag'',''hv1Value'');';
UI(5).Statement = 'set(hv1Value,''String'',num2str(hv,''%.1f''));';


% Save all the structures to a .mat file. & auto start
filename = ['MatFiles\UFD_cgrp'];
save(filename);

% return
VSX

return
% 
 % External function definitions
%#EF-1
Dop = computeSVD(RData)
    global na;
    IQ = squeeze(RData);
    IQ = reshape(IQ, size(IQ, 1), size(IQ, 2), na, []);
    IQ = squeeze(mean(IQ, 3));
  
    ncut = 30;
    % SVD PowerDoppler
    IQ_signal = IQ;%(:,:,1:100);
    [nz, nx, nt] = size(IQ_signal);
    IQ_signal = reshape(IQ_signal, [nz*nx, nt]);
    cov_matrix = IQ_signal'*IQ_signal;
    [Eig_vect, Eig_val]= eig(cov_matrix);
    Eig_vect=fliplr(Eig_vect);
    Eig_val=rot90(Eig_val,2);
    M_A = IQ_signal*Eig_vect;
    skipped_eig_val =[1:ncut];
    IQF_tissu = M_A(:,skipped_eig_val)*Eig_vect(:,skipped_eig_val)';
    IQF_tissu = reshape(IQF_tissu, [nz, nx, nt]);
    IQ_signal = reshape(IQ_signal, [nz, nx, nt]);
    IQF_corrected = IQ_signal-IQF_tissu;
    
    Dop = mean(abs(IQF_corrected(:,:,:)).^2,3);
    Dop = (Dop - min(Dop(:)))./range(Dop(:));
    
    figure(56)
    imagesc(interp2(sqrt(Dop),2)), colormap hot(256)
    caxis([0.05 0.95])
return
%#EF-1


%#EF-2
saveIQData(varargin)
    persistent bloc_count;
    persistent dir_save;
    global na;
    if isempty(bloc_count)
        bloc_count = 1;
    end
    IQData = complex(varargin{1}, varargin{2});
    tmp = squeeze(IQData);
    tmp = reshape(tmp, size(tmp, 1), size(tmp, 2), na, []);
    tmp = squeeze(mean(tmp, 3));   % do the compounding (less disk space needed)
    IQ_re = real(tmp);
    IQ_im = imag(tmp);
    IQ = [IQ_re IQ_im];
    path_save = 'D:\'; %Data save path
    if bloc_count == 1
        dir_save = ['Loop_fUS_' datestr(now, 'yyyy-mm-dd@HH-MM-SS') '\'];
        mkdir([path_save dir_save])
    end
    file_name = sprintf('fUS_block_%.3d.bin', bloc_count);
    
    fid = fopen([path_save dir_save file_name],'w');   % fast save
    fwrite(fid,IQ, 'double');
    fclose(fid);
    bloc_count = bloc_count+1;
    %save([path_save dir_save file_name], 'IQ', '-v6');
%#EF-2

%#EF-3
dispBlockTime()
    persistent time_bloc;
    if isempty(time_bloc)
        time_bloc = tic;
    end
    T = toc(time_bloc) % absolute time (ms)
    D = datestr(now, 'HH-MM-SS.FFF') % fUS image interval (s)
    
    fidresult = fopen('D:\timeLog.txt', 'a'); % time log save path
    fprintf(fidresult,'%s  %g\r\n', D, T);
    fclose(fidresult);
%#EF-3
