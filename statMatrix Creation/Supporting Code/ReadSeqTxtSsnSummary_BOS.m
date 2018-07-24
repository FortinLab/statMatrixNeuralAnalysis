function [txtDta] = ReadSeqTxtSsnSummary_BOS

%%
gettingFiles = 1;
txtFiles = [];
fNum = 1;
while gettingFiles
    [tempFname, tempPath] = uigetfile('.txt');
    if tempFname == 0
        gettingFiles = 0;
    else
        txtFiles{fNum} = [tempPath tempFname];
        fNum = fNum+1;
    end
end

%%
ratName = cell(size(txtFiles));
sniffTime = nan(size(txtFiles));
bufferDur = nan(size(txtFiles));
trlData = cell(size(txtFiles));

for file = 1:length(txtFiles)
    fID = fopen(txtFiles{file});
    txt = 1;
    while txt~=-1
        txt = fgetl(fID);
        if txt==-1            
            trlData{file} = struct('TrialNum', trialNum(1:trl-1), 'SeqNum', seqNum(1:trl-1),...
                'Odor', odors(1:trl-1), 'Position', positions(1:trl-1),...
                'TransDist', transDist(1:trl-1),...
                'PokeDur', pokeDur(1:trl-1), 'PokeIn', pokeInTime(1:trl-1),...
                'PokeOut', pokeOutTime(1:trl-1), 'Perf', perf(1:trl-1),...
                'NumPokes', numPokes(1:trl-1),...
                'RatName', repmat(ratName(file), [1,trl-1]),...
                'TargDur', repmat({sniffTime(file)}, [1,trl-1]),...
                'Buffer', repmat({bufferDur(file)}, [1,trl-1]));
        else
            if isempty(txt)
                txt = '1';
            end
            if ~isempty(strfind(txt, 'RatName'))
                [~, tempEnd] = regexp(txt, 'RatName = ');
                ratName{file} = txt(tempEnd+1:end);
            elseif ~isempty(regexp(txt, '^SniffTime =', 'once'))
                sniffTime(file) = str2double(cell2mat(regexp(txt, '([0-9]*)\.([0-9]*)', 'match')));
            elseif ~isempty(strfind(txt, 'NumOfTrials'))
                numTrials = str2double(cell2mat(regexp(txt, '([0-9]*)', 'match')));
                trialNum = cell(1,numTrials);
                seqNum = cell(1,numTrials);
                odors = cell(1,numTrials);
                positions = cell(1,numTrials);
                transDist = cell(1,numTrials);
                pokeDur = cell(1,numTrials);
                pokeInTime = cell(1,numTrials);
                pokeOutTime = cell(1,numTrials);
                perf = cell(1,numTrials);
                numPokes = repmat({1}, [1,numTrials]);
                trl = 1;
                pos = 1;
                seq = 1;
            elseif ~isempty(strfind(txt, 'OutInBuffer'))
                bufferDur(file) = str2double(cell2mat(regexp(txt, '([0-9]*)', 'match')));
            elseif ~isempty(regexp(txt, 'Odor ([0-9]*) delivered', 'once'))
                odors{trl} = str2double(cell2mat(regexp(txt, ' ([0-9]*) ', 'match')));
                positions{trl} = pos;
                transDist{trl} = pos - odors{trl};
                seqNum{trl} = seq;
            elseif ~isempty(regexp(txt, '^TimeInPort', 'once')) || ~isempty(regexp(txt, '^Time in port', 'once'))
                pokeDur{trl} = str2double(regexp(txt, '([0-9]*)\.([0-9]*)', 'match', 'once'));
                if ~isempty(strfind(txt, 'REWARD')) || ~isempty(strfind(txt, 'CORRECT'))
                    perf{trl} = 1;
                elseif ~isempty(strfind(txt, 'No reward'))
                    perf{trl} = 0;
                else
                    error('No performance indicator on poke dur line... check logic');
                end
            elseif ~isempty(strfind(txt, 'One poke ignored'))
                numPokes{trl} = numPokes{trl}+1;
            elseif ~isempty(regexp(txt, 'Rat came out ([0-9]*) times during odor presentation', 'once'))
                numPokes{trl} = numPokes{trl} + str2double(cell2mat(regexp(txt, '([0-9]*)', 'match')));
            elseif ~isempty(strfind(txt, 'Rat was in at'))
                pokeInTime{trl} = str2double(cell2mat(regexp(txt, '([0-9]*)\.([0-9]*)', 'match')));
            elseif ~isempty(strfind(txt, 'Rat was out at'))
                pokeOutTime{trl} = str2double(cell2mat(regexp(txt, '([0-9]*)\.([0-9]*)', 'match')));
            elseif ~isempty(regexp(txt, 'Trial ([0-9]*) completed', 'once'))
                trialNum{trl} = str2double(cell2mat(regexp(txt, '([0-9]*)', 'match')));
                trl = trl+1;
                pos = pos+1;
            elseif ~isempty(strfind(txt, 'Sequence completed correctly.'))
                pos = 1;
                seq = seq+1;
            elseif ~isempty(strfind(txt, 'Error in the sequence')) || ~isempty(regexp(txt, 'Sequence presented ([0-9]*) times', 'once'))
                pos = 1;
                seq = seq+1;
            end
        end
    end
    fclose(fID);
end

%%
txtDta = trlData{1};
for fl = 2:length(trlData)
    curSsnDta = trlData{fl};
    for trl = 1:length(curSsnDta)
       curSsnDta(trl).TrialNum = curSsnDta(trl).TrialNum + txtDta(end).TrialNum;
       curSsnDta(trl).SeqNum = curSsnDta(trl).SeqNum + txtDta(end).SeqNum;
    end
    txtDta = [txtDta, curSsnDta]; %#ok<AGROW>
end

%%
uisave('txtDta', '_txtDta');
            