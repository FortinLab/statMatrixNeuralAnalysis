function [behavData] = CombineTrialLogVects_SM(behav1, behav2)
%% CombineTrialLogVects_SM
%   Function to combine logical vectors created in OrganizeTrialData_SM
%   using a window of [0 0] (to extract logical indices with timestamp
%   data). Main use will be to create logical vectors scaled to match the
%   trial size, i.e. behav1=pokeIn aligned, behav2=pokeOut aligned.
%
%% Check inputs
if ~(length(behav1)==length(behav2))
    error('Trial vectors are of unequal size, check inputs');
elseif sum([behav1.SequenceNum]-[behav2.SequenceNum])==0 ||...
        sum([behav1.Odor]-[behav2.Odor])==0 ||...
        sum([behav1.Position]-[behav2.Position])==0 ||...
        sum([behav1.PokeDuration]-[behav2.PokeDuration])==0 ||...
        sum([behav1.Performance]-[behav2.Performance])==0 ||...
        sum([behav1.TranspositionDistance]-[behav2.TranspositionDistance])==0 ||...
        nansum([behav1.ItemItemDistance]-[behav2.ItemItemDistance])==0
    behavData = behav1;
    for trl = 1:length(behavData)
        b1Log = behav1(trl).TrialLogVect;
        b2Log = behav2(trl).TrialLogVect;
        behavData(trl).TrialLogVect(find(b1Log,1,'first'):find(b2Log,1,'last')) = true;
    end
end