%% Runtime variables
piWindow = [-0.5 0.5];
poWindow = [-0.5 0.5];
slideWin = 200;
dsRate = 5;
phaseBins = -pi:pi/3.5:pi;
cLim = [0 0.01];
spectRanges = [[4,12];...
    [16, 32];...
    [40 59];...
    [61 100]];

piWindowTS = piWindow(1):(dsRate*mode(diff(tsVect))):piWindow(2);
poWindowTS = poWindow(1):(dsRate*mode(diff(tsVect))):poWindow(2);

%% Create the spikeMtx: the multiunit signal from every tetrode by aggregating the cut and uncut units together.
tempTSvect = tsVect;
tempTSvect(end+1) = tsVect(end)+(mode(diff(tsVect)));
spikeMtx = nan(length(tsVect), length(osmNeural));
lfpInstPhase = nan(length(tsVect), length(osmNeural), size(spectRanges,1));

for t = 1:length(osmNeural)
    tempTetSpks = osmNeural(t).Spikes;
    tempSpikes = [];
    if ~isempty(tempTetSpks)
        uniNms = fieldnames(tempTetSpks);
        if length(uniNms)~=1
            tempTetLFP = osmNeural(t).LFP;
            for f = 1:size(spectRanges,1)
                lfpInstPhase(:,t,f) = SimpleHilbFilt(tempTetLFP, 1/mode(diff(tsVect)), spectRanges(f,:));
            end
            for u = 2:length(uniNms)
                tempSpikes = [tempSpikes; tempTetSpks.(uniNms{u})]; %#ok<AGROW>
            end
            spikeMtx(:,t) = histcounts(tempSpikes, tempTSvect);
        end
    end
end
tetNames = {osmNeural.TetName};
tetNames(isnan(spikeMtx(1,:))) = [];
lfpInstPhase(:,isnan(spikeMtx(1,:)),:) = [];
for ph = 1:size(lfpInstPhase,3)
    lfpInstPhase(:,:,ph) = [lfpInstPhase(:,end,ph), lfpInstPhase(:,1:end-1, ph)];       % Rotate LFP to other tetrodes to break LFP/spiking relationship
%     lfpInstPhase(:,:,ph) = repmat(lfpInstPhase(:,1,ph), [1,size(lfpInstPhase,2)]);      % Or... use a common tetrode's LFP signal
end
spikeMtx(:,isnan(spikeMtx(1,:))) = [];

%% Bin spiking on a per trial basis
piSpikeMtx = nan(length(piWindowTS)-1, length(tetNames), length(osmBehavior));
poSpikeMtx = nan(length(poWindowTS)-1, length(tetNames), length(osmBehavior));
trlSpikeMtx = nan(length(piWindowTS) + length(poWindowTS) -2, length(tetNames), length(osmBehavior));

piSpikeMtxPhase = repmat({piSpikeMtx}, [size(spectRanges,1), length(phaseBins)-1]);
poSpikeMtxPhase = repmat({piSpikeMtx}, [size(spectRanges,1), length(phaseBins)-1]);
trlSpikeMtxPhase = repmat({[piSpikeMtx;poSpikeMtx]}, [size(spectRanges,1), length(phaseBins)-1]);

lfpInstPhaseTrlMtx = repmat({nan(length(piWindowTS) + length(poWindowTS) -2, length(tetNames), length(osmBehavior))}, [1 size(spectRanges,1)]);

for trl = 1:length(osmBehavior)
    pokeInWindowLog = tsVect>osmBehavior(trl).PortEntryTime + piWindow(1) - (slideWin/2*mode(diff(tsVect))) & tsVect<=osmBehavior(trl).PortEntryTime + piWindow(2) + (slideWin/2*mode(diff(tsVect)));
    pokeOutWindowLog = tsVect>osmBehavior(trl).PortExitTime + poWindow(1) - (slideWin/2*mode(diff(tsVect))) & tsVect<=osmBehavior(trl).PortExitTime + poWindow(2) + (slideWin/2*mode(diff(tsVect)));    
    
    for tet = 1:length(tetNames)
        curPIspikes = spikeMtx(pokeInWindowLog,tet);
        tempSpksPI = conv(curPIspikes, ones(1,slideWin)./(slideWin/1000), 'same');
        piSpikeMtx(:,tet,trl) = downsample(tempSpksPI((slideWin/2)+1:end-(slideWin/2)),dsRate);
        
        curPOspikes = spikeMtx(pokeOutWindowLog,tet);
        tempSpksPO = conv(curPOspikes, ones(1,slideWin)./(slideWin/1000), 'same');
        poSpikeMtx(:,tet,trl) = downsample(tempSpksPO((slideWin/2)+1:end-(slideWin/2)),dsRate);
        
        trlSpikeMtx(:,tet,trl) = [piSpikeMtx(:,tet,trl);poSpikeMtx(:,tet,trl)];
        for freq = 1:size(spectRanges,1)
            curTetFreqHilbPI = lfpInstPhase(pokeInWindowLog,tet,freq);
            curTetFreqHilbPO = lfpInstPhase(pokeOutWindowLog,tet,freq);
            lfpInstPhaseTrlMtx{freq}(:,tet,trl) = [downsample(curTetFreqHilbPI((slideWin/2)+1:end-(slideWin/2)), dsRate);...
                downsample(curTetFreqHilbPO((slideWin/2)+1:end-(slideWin/2)), dsRate)];
            for phase = 2:length(phaseBins)
                tempPIspks = curPIspikes;
                curPhaseLogPI = curTetFreqHilbPI>=phaseBins(phase-1) & curTetFreqHilbPI<phaseBins(phase);
                tempPIspks(~curPhaseLogPI) = 0;
                tempPIspks = conv(tempPIspks, ones(1,slideWin)./(slideWin/1000), 'same');
                piSpikeMtxPhase{freq,phase-1}(:,tet,trl) = downsample(tempPIspks((slideWin/2)+1:end-(slideWin/2)),dsRate);
                                
                tempPOspks = curPOspikes;
                curPhaseLogPO = curTetFreqHilbPO>=phaseBins(phase-1) & curTetFreqHilbPO<phaseBins(phase);
                tempPOspks(~curPhaseLogPI) = 0;
                tempPOspks = conv(tempPOspks, ones(1,slideWin)./(slideWin/1000), 'same');
                poSpikeMtxPhase{freq,phase-1}(:,tet,trl) = downsample(tempPOspks((slideWin/2)+1:end-(slideWin/2)),dsRate);
                
                trlSpikeMtxPhase{freq,phase-1}(:,tet,trl) = [piSpikeMtxPhase{freq,phase-1}(:,tet,trl);poSpikeMtxPhase{freq,phase-1}(:,tet,trl)];
            end
        end                
    end
    fprintf('Trial %i\n', trl)
end

%% Now arrange the logicals & likelihoods
% Use these to select trials
trlPos = [osmBehavior.TrialPosition];
trlOdr = [osmBehavior.TrialOdor];
isTrlLog = (trlPos - trlOdr)==0;
perfTrlLog = [osmBehavior.PerformanceLog]==1;

fisLog = false(1,length(osmBehavior));
for t = 4:length(osmBehavior)
    fisLog(t) = sum([osmBehavior(t-3:t).TrialPosition] == 1:4)==4 &&...
        sum([osmBehavior(t-3:t).TrialOdor] == 1:4)==4 &&...
        sum([osmBehavior(t-3:t).PerformanceLog])==4;
    if fisLog(t) == true
        fisLog(t-3:t) = true;
    end
end

% Arrange the beast of all likelihood matrices...
freqBandIDvect = repmat({nan(length(piWindowTS) + length(poWindowTS) -2,1)}, [size(spectRanges,1), length(phaseBins)-1]);
phaseBinIDvect = repmat({nan(length(piWindowTS) + length(poWindowTS) -2, 1)}, [size(spectRanges,1), length(phaseBins)-1]);
for f = 1:size(spectRanges,1)
    for p = 1:length(phaseBins)-1
        freqBandIDvect{f,p} = ones(size(freqBandIDvect{f,p}))*f;
        phaseBinIDvect{f,p} = ones(size(phaseBinIDvect{f,p}))*p;
    end
end
freqBandIDvect = repmat(cell2mat(freqBandIDvect(:)), [4,1]);
phaseBinIDvect = repmat(cell2mat(phaseBinIDvect(:)), [4,1]);

odrIDphaseVect = nan((length(piWindowTS) + length(poWindowTS) - 2) * size(spectRanges,1) * (length(phaseBins)-1)*4,1);
odrIDvect = zeros(1,(length(piWindowTS)-1 + length(poWindowTS)-1)*4);
for o = 1:4
    odrIDvect((o-1)*(length(piWindowTS)+length(poWindowTS)-2)+1:(o*(length(piWindowTS)+length(poWindowTS)-2))) = o;
    odrIDphaseVect((o-1)*((length(piWindowTS) + length(poWindowTS) - 2) * size(spectRanges,1) * (length(phaseBins)-1))+1:...
        (o*((length(piWindowTS) + length(poWindowTS) - 2) * size(spectRanges,1) * (length(phaseBins)-1)))) = o;
end

trlSpikePhaseVect = cell2mat(trlSpikeMtxPhase(:));

%% Extract FIS Trials to create likelihood matrix
fisSpikeMtx = nan(length(odrIDvect), length(tetNames), sum(fisLog)/4);
fisSpikePhaseMtx = nan(length(odrIDphaseVect), length(tetNames), sum(fisLog)/4);
for odr = 1:4
    tempFISodrVect = find(fisLog & trlOdr==odr);
    for tet = 1:length(tetNames)
        for trl = 1:length(tempFISodrVect)
            fisSpikeMtx(odrIDvect==odr,tet,trl) = trlSpikeMtx(:,tet,tempFISodrVect(trl));
            fisSpikePhaseMtx(odrIDphaseVect==odr,tet,trl) = trlSpikePhaseVect(:,tet,tempFISodrVect(trl));
        end
    end
end

%% Decode stuff
fisPost = nan(size(fisSpikeMtx,1), size(fisSpikePhaseMtx,1), size(fisSpikeMtx,3));
for seq = 1:size(fisSpikeMtx,3)
    curObsv = fisSpikeMtx(:,:,seq);
    tempLikely = fisSpikePhaseMtx;
    tempLikely(:,:,seq) = [];
    fisPost(:,:,seq) = CalcStaticBayesPost(mean(tempLikely,3), curObsv, slideWin);
end

figure;
fisFREQdecode = DecodeBayesPost(fisPost, freqBandIDvect);
subplot(3,1,1)
for f = 1:size(spectRanges,1)
    plot(mean(fisFREQdecode==f,2));
    hold on;
end
fisPHASEdecode = DecodeBayesPost(fisPost, phaseBinIDvect);
subplot(3,1,2)
for p = 1:length(phaseBins)-1
    plot(mean(fisPHASEdecode==p,2));
    hold on;
end
fisODRdecode = DecodeBayesPost(fisPost, odrIDphaseVect);
subplot(3,1,3)
for o = 1:4
    plot(mean(fisODRdecode==o,2));
    hold on;
end
    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functionarinos!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SimpleHilbFilt
function zeroTroughPhase = SimpleHilbFilt(data, samp, window)
    Wn_FRange = [window(1)/(samp/2) window(2)/(samp/2)]; % normalized by the nyquist frequency
    [bFRange, aFRange] = butter(3,Wn_FRange);
    filtered = filtfilt(bFRange,aFRange,data);

    zeroTroughPhase = atan2(imag(hilbert(filtered*-1)), real(hilbert(filtered*-1)));
end


%% CalcStaticBayesPost
function post = CalcStaticBayesPost(likely, obsv, binSize)
post = nan(size(obsv,1), size(likely,1), size(obsv,3));
for trl = 1:size(obsv,3)
    for t = 1:size(obsv,1)
        p = nan(size(likely));
        curPopVect = obsv(t,:,trl)*(binSize/1000);
        curPopFact = factorial(curPopVect);
        for u = 1:size(likely,2)
            curAvgUniFR = likely(:,u);
            p(:,u) = (((binSize/1000).*curAvgUniFR).^curPopVect(u))./curPopFact(u);
        end        
        pp = prod(p,2);
        ee = exp(-((binSize/1000)*sum(likely,2)));
        tempPost = pp.*ee;
        post(t,:,trl) = tempPost./sum(tempPost);
    end
end
end
    
%% DecodeBayesPost
function decode = DecodeBayesPost(post, id)
% Assumes post is in the structure of ObservTime X LikelyTime X Trial
decode = nan(size(post,1),size(post,3));
for o = 1:size(post,3)
    for d = 1:size(post,1)
        if ~isnan(post(d,1,o))
            decode(d,o) = id(find(post(d,:,o)==max(post(d,:,o)),1,'first'));
        end
    end
end
end