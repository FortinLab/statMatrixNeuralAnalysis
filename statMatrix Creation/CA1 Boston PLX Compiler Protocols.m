%% Step 1: Create the txtData file
[txtDta] = ReadSeqTxtSsnSummary_BOS;

%% Step 2: 
% Navigate to the appropriate file directory
% And identify the plxFile being used. Make sure that if the session's data
% is composed of multiple plx files you merge them appropriately using the
% PlexUtil program.

% plxFile = 'SuperChris-2-12-09.plx'
plxFile = 'Stella-2-12-2009.plx'
% plxFile = 'Barat-11-06-2008Skips_mrg.plx'
% plxFile = 'Buchanan4-20-withskips_mrg.plx'
% plxFile = 'Mitt_July18_5odorswithSkips.plx'

% Run the plxData compiler.
%   -It will generate a text file with the suffix "_MergeSummary_####" 
%       NOTE: the numbers are arbitrary but preserved across runs of the
%       same data constructed in the SummarizePLXabbr_BOS2 function. So all
%       the txt file outputs will have the same number.
%
% IF the order of the odors in the text files and plx files do not match
%   a list of the odors in both files will be spit out in the text file.
%   
%   -Review the resulting text file and rename it with the
%       suffix "_OdorListMismatch.txt" so that the file trail is preserved.
%   -If the lists don't match, you will have to manually resolve the mismatch
%       so that the text files match the plx file. 
%       
%       -BE CAREFUL about removing the text file sessions. 
%       -MOST of the discrepancies occur towards the end of the text files 
%           where acquisition on the plexon rig was halted but the
%           matlab code continued to run. 
%       -This can lead to ambiguity, as to which trials should be kept
%           if there are multiple sessions run for an animal during the day
%       -WHEN IN DOUBT remove trials at the END of a session to account for 
%           missing trials.
% 
[plxData] = SummarizePLXabbr_BOS3(plxFile);
% The output of this function creates a figure that's saved as
%   "[filename]_PokeDur_Summary.fig"
% When the order of odor presentations is resolved, the output of the
%   "MergeSummary" file will contain a list of times there were "errors" in
%   the plexon code. MOST of these errors will stem from buffer related
%   issues, however some will potentially change how the trial SHOULD have
%   been interpreted. You will manually curate them in the next section.
if isempty(plxData)
    return
end
fclose all;
% These functions plot the poke data for the session as asterixes and lines
% to show the poking behavior:
%   For All Trial: (Save this figure with the suffix "_TrialPokes_ALL")
PlotErrorTrialsPLXtxtMerge(plxData)
%   For only trials flagged as "errors": "Save this figure with the suffix
%   "_TrialPokes_Errors")
PlotAllTrialsPLXtxtMerge(plxData)

%% Step 3: Curate the plxData file
% Here is where the "MergeSummary" output is put to work. This function
% also creates a text file output "_CuratePLX_####". NOTE: as mentioned
% above, the number will be preserved across runs to maintain
% continuity through a processing instance.
%
% This function will step through the "error" trials and make you decide on
% each trial whether to use the Plexon timestamp, the text file timestamp,
% OR remove the trial entirely.
%   Reviewing these trials should be done by:
%       A: Consulting the "MergeSummary" text file as that contains
%           information about the error
%       B: Reviewing the video file using the CinePlexMarkup program. The
%           "MergeSummary" file has the poke initiation time listed in it
%           for each of the "Error" trials.
%           NOTE:
%           with CinePlexMarkup you will need to have both the .plx file(s)
%           as well as the .avi file(s) for the session on hand. 
%           CinePlexMarkup will then combine them into a .cpj file. 
%           ALSO NOTE: 
%           If there are multiple .plx/.avi files for a session you will
%           (likely) need to adjust the poke in timestamps if they are
%           coming from multiple files (i.e. if you're looking for a trial
%           at 1350 and the session had two files, one that ended at 1300,
%           you will then need to look at the second file from that day at
%           timestamp 50 (1350 - 1300)... **I actually haven't had to do
%           this yet, the subtraction may not be an accurate approach if
%           the merging of the .plx files modified the timestamps**
[plxData] = CuratePLXssnDataTxtMerge(plxData);
% This function save the resulting plxData into a separate file. USE THAT 
%   FILE FOR CREATING THE STAT MATRIX (below)
fclose all;
% This file will also output a figure showing the NEW timestamps following
% curation. Trials you removed will be missing (if that's not the case
% review the CuratePLX file to make sure you did it right). If there's an
% error, delete the CuratePLX as well as the _plxData files and run this
% section again.

%% Step 4: Create the statMatrix files
% Still need to tweak this statMatrix creator function. It's all commented
% but I need to add in the txt file output for this to preserve the paper
% trail. I also may modify the order of the dialogue boxes so they make
% more sense.

StatMatrixCreator