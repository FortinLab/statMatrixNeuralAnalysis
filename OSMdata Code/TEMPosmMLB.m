%% Runtime variables
piWindow = [-0.5 0.5];
poWindow = [-0.5 0.5];
binSize = 200;
dsRate = 5;
cLim = [0 0.01];

piWindowTS = piWindow(1):(dsRate*mode(diff(tsVect))):piWindow(2);
poWindowTS = poWindow(1):(dsRate*mode(diff(tsVect))):poWindow(2);

%% Create the spikeMtx: the multiunit signal from every tetrode by aggregating the cut and uncut units together.
tempTSvect = tsVect;
tempTSvect(end+1) = tsVect(end)+(mode(diff(tsVect)));
spikeMtx = nan(length(tsVect), length(osmNeural));
for t = 1:length(osmNeural)
    tempTetSpks = osmNeural(t).Spikes;
    tempSpikes = [];
    if ~isempty(tempTetSpks)
        uniNms = fieldnames(tempTetSpks);
        for u = 1:length(uniNms)
            tempSpikes = [tempSpikes; tempTetSpks.(uniNms{u})]; %#ok<AGROW>
        end
        spikeMtx(:,t) = histcounts(tempSpikes, tempTSvect);
    end
end
tetNames = {osmNeural.TetName};
tetNames(isnan(spikeMtx(1,:))) = [];
spikeMtx(:,isnan(spikeMtx(1,:))) = [];

%% Bin spiking on a per trial basis
piSpikeMtx = nan(length(piWindowTS)-1, length(tetNames), length(osmBehavior));
poSpikeMtx = nan(length(poWindowTS)-1, length(tetNames), length(osmBehavior));
trlSpikeMtx = nan(length(piWindowTS) + length(poWindowTS) -2, length(tetNames), length(osmBehavior));
for trl = 1:length(osmBehavior)
    pokeInWindowLog = tsVect>osmBehavior(trl).PortEntryTime + piWindow(1) - (binSize/2*mode(diff(tsVect))) & tsVect<=osmBehavior(trl).PortEntryTime + piWindow(2) + (binSize/2*mode(diff(tsVect)));
    pokeOutWindowLog = tsVect>osmBehavior(trl).PortExitTime + poWindow(1) - (binSize/2*mode(diff(tsVect))) & tsVect<=osmBehavior(trl).PortExitTime + poWindow(2) + (binSize/2*mode(diff(tsVect)));    
    
    for tet = 1:length(tetNames)
        tempSpksPI = conv(spikeMtx(pokeInWindowLog,tet), ones(1,binSize)./(binSize/1000), 'same');
        tempSpksPI = tempSpksPI((binSize/2)+1:end-(binSize/2));
        piSpikeMtx(:,tet,trl) = downsample(tempSpksPI,dsRate);
        
        tempSpksPO = conv(spikeMtx(pokeOutWindowLog,tet), ones(1,binSize)./(binSize/1000), 'same');
        tempSpksPO = tempSpksPO((binSize/2)+1:end-(binSize/2));
        poSpikeMtx(:,tet,trl) = downsample(tempSpksPO,dsRate);
        
        trlSpikeMtx(:,tet,trl) = [piSpikeMtx(:,tet,trl);poSpikeMtx(:,tet,trl)];
    end
    fprintf('Trial %i\n', trl)
end

%% Create logical masks
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

% Use these to navigate the FIS trial spike matrix (created below)
piIDvect = repmat([true(1,length(piWindowTS)-1), false(1,length(poWindowTS)-1)], [1,4]);
poIDvect = repmat([false(1,length(piWindowTS)-1), true(1,length(poWindowTS)-1)], [1,4]);
odrIDvect = zeros(1,(length(piWindowTS)-1 + length(poWindowTS)-1)*4);
for o = 1:4
    odrIDvect((o-1)*(length(piWindowTS)+length(poWindowTS)-2)+1:(o*(length(piWindowTS)+length(poWindowTS)-2))) = o;
end
timIDvect = repmat(1:length(piWindowTS)+length(poWindowTS)-2, [1,4]);

%% Extract FIS Trials to create likelihood matrix
fisSpikeMtx = nan(length(odrIDvect), length(tetNames), sum(fisLog)/4);
for odr = 1:4
    tempFISodrVect = find(fisLog & trlOdr==odr);
    for tet = 1:length(tetNames)
        for trl = 1:length(tempFISodrVect)
            fisSpikeMtx(odrIDvect==odr,tet,trl) = trlSpikeMtx(:,tet,tempFISodrVect(trl));
        end
    end
end

%% Perform leave 1 out decoding on FIS data
fisPost = nan(size(fisSpikeMtx,1), size(fisSpikeMtx,1), size(fisSpikeMtx,3));
odrDecode = nan(size(fisSpikeMtx,1),size(fisSpikeMtx,3));
for seq = 1:size(fisSpikeMtx,3)
    curObsv = fisSpikeMtx(:,:,seq);
    tempLikely = fisSpikeMtx;
    tempLikely(:,:,seq) = [];
    fisPost(:,:,seq) = CalcStaticBayesPost(mean(tempLikely,3), curObsv, binSize);
end
fisODRdecode = DecodeBayesPost(fisPost, odrIDvect);
fisODRdecoded = nan(4,size(fisODRdecode,1));
for o = 1:4
    fisODRdecoded(o,:) = smooth(mean(fisODRdecode == o,2),10);
end
fisLATdecode = DecodeBayesPost(fisPost, timIDvect);
for seq = 1:size(fisLATdecode,2)
    fisLATdecode(:,seq) = fisLATdecode(:,seq) - timIDvect';
end
lagMean = smooth(mean(fisLATdecode,2),10);
lagVar = std(fisLATdecode,0,2);


%% Plot the MLB summary
figure; 
% Plot the probability density
subplot(5,5,[1:4, 6:9, 11:14, 16:19]);
imagesc(mean(fisPost,3), [0 0.01]);
set(gca, 'ydir', 'normal');
colormap jet;
% Plot odor decoding through the sequence
subplot(5,5,21:24);
plot(fisODRdecoded(1,:), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(fisODRdecoded(2,:), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(fisODRdecoded(3,:), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(fisODRdecoded(4,:), 'color', [128/255 66/255 151/255], 'linewidth', 1);
% Plot temporal decoding through the sequence
subplot(5,5,5:5:20)
plot(lagMean, 1:size(fisLATdecode,1), 'k');
hold on;
plot([0 0], get(gca, 'ylim'), ':k');
patch('XData', [lagMean+lagVar; flipud(lagMean-lagVar)],...
    'YData', [1:length(lagMean), length(lagMean):-1:1], 'FaceColor', 'k', 'FaceAlpha', .3, 'edgecolor', 'none');


%% Functions %%

%%
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
    
%%
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