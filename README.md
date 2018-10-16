# statMatrixNeuralAnalysis
**************************************************************
# Overview of the statMatrix
**************************************************************
The statMatrix data structure is a standard data matrix structure used by the Fortin Lab. It was designed to provide a standard organization for both neural and behavioral data derived from experimental sessions. Each file contains two workspace variables, 1) the main data structure and 2) column identifiers for the main data structure. The main data structure consists of an NxM double matrix consisting of N rows, determined by the sample rate of the LFP signal for the recording session (standard = 1kHz), and M columnns corresponding to different information sources. For each recording session a separate statMatrix file is produced for each recording source, i.e. tetrode as well as a separate file with containing the behavior data from the session.

For data stored from each recording source, i.e. tetrode, the workspace variables saved in those files are the 'statMatrix' (main data matrix) and the 'statMatrixColIDs' (column identifiers), whereas behavioral data contains the 'behavMatrix' (main data matrix) and the 'behavMatrixColIDs' (column identifiers).

The statMatrix is organized with rows indexed to the LFP sampleRate. This makes it easy to associate spiking activity and behavioral events to LFP signals with minimal loss of precision. As the LFP data is either collected directly at 1kHz/s or downsampled to that frequency, the loss of precision, i.e. associating a spike/event to one ms or another, is trivial, especially since most analysis is done using time aggregated spiking (spk/s).

***************************************************************
# Working with the statMatrix
***************************************************************
____________________________________________
### statMatrix Columns Organization
____________________________________________
**Sequence Task**
* **Timebin:** Timestamps pulled from the LFP trace
* **LFP Data:** Multiple columns consisting of the Raw LFP trace as well as bandpass filtered traces for band-specific analysis. Each frequency range contains two columns, one indicating the voltage value for that trace e.g. "_RAW" or "_Theta" as well as a column of phase values appended with *"_HilbVals," e.g. "_RAW_HilbVals"* or *"_Theta_HilbVals."* 
* **Unit Data:** Logical vector (there shouldn't be any 2s...) indicating individual unit spiking activity. 1s indicate if the unit spiked during that time bin. 
____________________________________________
### behavMatrix Columns Organization
____________________________________________
**Sequence Task**
* **'Odor\[1-X]'**: Columns with logical 1 indicating when odor was delivered (no flag or indicator when odor presentation was terminated, assume it was at port withdrawal or trial feedback, whichever came first)
* **'Position\[1-X]'**: Columns with logical 1 indicating what sequence position the trial occured during. Indexed with odor delivery.
* **'InSeqLog'**: Column with values of \[1,0,-1]; 1 = InSeq trial, 0 = Nothing/Filler, -1 = OutSeq trial. Use inSeq==1 and outSeq==-1 to create logical trial vectors; use 'trials==(abs(InSeqLogColumn#)==1)' or something like it can be used to pull out trial start indices. Indexed with odor delivery.
* **'PerformanceLog'**: Column with values \[1,0,-1]; 1 = Correct trial, 0 = Nothing/Filler, -1 = Incorrect trial. Can be used like 'InSeqLog' to pull out correct/incorrect trials and/or trial start indices. Indexed with odor delivery.
* **'PokeEvents'**: Column with values \[1,0,-1]; 1 = Port Entry, 0 = Nothing/Filler, -1 = Port Withdrawal. To identify port entry relative to trial start identify the last port entry (PokeEvents==1) prior to odor delivery; likewise to identify port withdrawal identify the first port withdrawal (PokeEvents==-1) following odor delivery.
* **'FrontReward'**: Column with logical 1 indicating when reward was given at the front of the maze.
* **'BackReward'**: Column with logical 1 indiciating when reward was given at the back of the maze.
* **'XvalRatMazePosition'**: Column indexing the rat's position within the maze along the long axis of the maze.
* **'YvalRatMazePosition'**: Column indexing the rat's position within the maze along the short axis of the maze. **NOTE** The motion capture system sample rate is lower than the time bins used to organize the statMatrix (~30Hz vs 1kHz), all non-zero positions are actual position values, position \[0,0] is out of the maze and the rat never went there.

****************************************************************
# List of statMatrix Functions
****************************************************************
NOTE Any modifications made to tailor code to processing a different file structure or data set should be saved as a new file and apppropriately named and commented to reflect that.

____________________________________________
### statMatrix Creation Functions
____________________________________________
The following functions create statMatrix data files for each tetrode recorded from during the training session. Note that not all tetrodes will have units but they will all have LFP data recorded from them. Currently these functions create the traditional per-tetrode statMatrix files (i.e. timestamp, LFP, Unit, Behavior sections) but they also output a separate file with the behavioral section of the statMatrix saved separately. This behavMatrix file is identical to the Behavior section of the statMatrix and is saved with a corresponding behavMatrixColIDs variable for column based indexing.

* **CreateStatMatrixFromPLX2.m**:
The initial statMatrix creation function. Designed to extract data from a plexon (.plx) session data file. Originally created for use with the .plx files recorded by NJF in Boston. The main difference between this code and other versions is the behavioral analysis code associated with it (SummarizePLXabbr_BOS.m). The rest of the code should be identical to other variants. NOTE: This code is currently legacy as the integration of MountainSort into the pre-processing stages necessitated creation of new compilation code.

* **CreateStatMatrixFromPLX2irvine.m**:
Variation of the original statMatrix creation function CreateStatMatrixFromPLX2.m designed to extract data from .plx session data files recorded at UCI. NOTE: This code is currently legacy following integration of MountainSort as a pre-processing stage.

* **CreateStatMatrixFromPLX_MS.m**:
Variation of the original statMatrix creation function designed to extract data from .plx files produced following MountainSort pre-processing. NOTE: This code is currently tailored for use with the Boston data (i.e. it's currently coded to work with SummarizePLXabbr_BOS). To convert it to working with UCI files it should be as simple as swapping out the Behavioral analysis function that needs to be tested before being done and will necessitate creation of a new file and validation before inclusion into the code set.

____________________________________________
### statMatrix Organization Functions
____________________________________________

* **EnsembleCompilation_SM.m**:
Code to extract all the unit data from all the individual per-tetrode statMatrix files and concatenate them together to create a single matrix. The data retains its original statMatrix organization in the variable called 'ensembleMatrix' and is paired with a column ID vector containing the corresponding unit names 'ensembleMatrixColIDs' for logical indexing to select or remove units.

* **OrganizeTrialData_SM.m**:
Code to organize the behavior data into a 1xN structure variable where each index corresponds to a session trial. Compiles information about each trial as subfields at each index and creates a logical vector for that trial period that can be used to extract neural data that occurred during that trial. **Needs to be commented**.

* **ExtractTrialData_SM.m**:
Code to use the trial period logical vector created in OrganizeTrialData_SM.m to extract neural data stored in matrices.

* **OrganizeStatMatrixByTrial.m**:
Obsolete! Feel free to review for example code. Will be removed... eventually.

************************************************************************
# List of Required Functions/Toolboxes
************************************************************************
* *Plexon Offline Files SDK*: 
Toolbox created by Plexon to analyze .plx files in Matlab. Link [here](https://plexon.com/wp-content/uploads/2017/08/OmniPlex-and-MAP-Offline-SDK-Bundle_0.zip)
* *Circular Statistics Toolbox*:
Toolbox from the matlab file exchange to perform circular (directional) statistics. Link [here](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-)
