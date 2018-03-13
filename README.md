# statMatrixNeuralAnalysis
**************************************************************
# Overview of the statMatrix
**************************************************************
The statMatrix data structure is a standard data matrix structure developed in the Fortin Lab. It was designed to contain the neural and behavioral data from an experimental session. Each file contains two workspace variables. The statMatrix data structure, an NxM double matrix consisting of N rows corresponding to each LFP time sample and M columnns corresponding to different information sources. The second variable is the statMatrixColIDs which is a 1xM cell vector containing strings indicating the identity of each column. The columns of the statMatrix is generally organized as follows:

Column Order 
(Note: index value may vary across files/experiments, numbers here reflect general order, not necessarily the corresponding index in the matrix):
1) Timebin: Timestamps pulled from the LFP trace
2) LFP Data: Multiple columns consisting of the Raw LFP trace as well as bandpass filtered traces for band-specific analysis. Each frequency range contains two columns, one indicating the voltage value for that trace e.g. "_RAW" or "_Theta" as well as a column of phase values appended with *"_HilbVals," e.g. "_RAW_HilbVals"* or *"_Theta_HilbVals."* (*Specifics for different experimental conditions should be added below*)
3) Unit Data: Logical column vector ^(there shouldn't be any 2s...)^ indicating individual unit spiking activity. 1s indicate if the unit spiked during that time bin. 
4) Behavior Data: Mostly logical vectors indicating when behavioral events occurred during the session. (*Specifics for different experimental conditions should be added below*)

***********************************************************
# Creating the statMatrix
***********************************************************
The statMatrix is made in Matlab using custom made m-files. Unique .m conversion files are used for constructing the statMatrix depending on the configuration used to acquire the behavioral and neurophysiological data. Though the neural data organization for .plx (Plexon) or .spikes/.continuous (OpenEphys) are standardized by file type, the behavioral timestamps associated with them (events channels for .plx and ADC .continuous channels) are not always so standardized. Therefore, to properly extract behavioral timestamps you need to be sure you are using the correct behavioral analysis code for the data based on the rig it was collected on. I will probably make an index of these associates at some point but for now there's too much to flesh out here.

As mentioned above the statMatrix is organized with rows indexed to the LFP sampleRate. This makes it easy to associate spiking activity and behavioral events to LFP signals with minimal loss of precision. As the LFP data is either collected directly at 1kHz/s or downsampled to that frequency, the loss of precision, i.e. associating a spike/event to one ms or another, is trivial, especially since most analysis is done on time aggregated measures (spk/s, spk per 


***************************************************************
# Working with the statMatrix
***************************************************************
Storage of data in statMatrix is advantageous... increased flexability... plug and play use... enables development of analysis/visualization tools that can be applied to any data set stored in that format...

Current status of statMatrix code... creation of per-tetrode data structure... creation of per-session data files

****************************************************************
# List of statMatrix Functions
****************************************************************
NOTE Any modifications made to tailor code to processing a different file structure or data set should be saved as a new file and apppropriately named and commented to reflect that.

____________________________________________
# statMatrix Creation Functions
____________________________________________
The following functions create statMatrix data files for each tetrode recorded from during the training session. Note that not all tetrodes will have units but they will all have LFP data recorded from them. Currently these functions create the traditional per-tetrode statMatrix files (i.e. timestamp, LFP, Unit, Behavior sections) but they also output a separate file with the behavioral section of the statMatrix saved separately. This behavMatrix file is identical to the Behavior section of the statMatrix and is saved with a corresponding behavMatrixColIDs variable for column based indexing.

---CreateStatMatrixFromPLX2.m---
The initial statMatrix creation function. Designed to extract data from a plexon (.plx) session data file. Originally created for use with the .plx files recorded by NJF in Boston. The main difference between this code and other versions is the behavioral analysis code associated with it (SummarizePLXabbr_BOS.m). The rest of the code should be identical to other variants. NOTE: This code is currently legacy as the integration of MountainSort into the pre-processing stages necessitated creation of new compilation code.

---CreateStatMatrixFromPLX2irvine.m---
Variation of the original statMatrix creation function CreateStatMatrixFromPLX2.m designed to extract data from .plx session data files recorded at UCI. NOTE: This code is currently legacy following integration of MountainSort as a pre-processing stage.

---CreateStatMatrixFromPLX_MS.m---
Variation of the original statMatrix creation function designed to extract data from .plx files produced following MountainSort pre-processing. NOTE: This code is currently tailored for use with the Boston data (i.e. it's currently coded to work with SummarizePLXabbr_BOS). To convert it to working with UCI files it should be as simple as swapping out the Behavioral analysis function that needs to be tested before being done and will necessitate creation of a new file and validation before inclusion into the code set.

____________________________________________
# statMatrix Organization Functions
____________________________________________
The following functions are written to work off the per-tetrode statMatrix data structures and organize them in different ways that may be useful for different analyses. 

---EnsembleCompilation_SM---
Code to extract all the unit data from all the individual per-tetrode statMatrix files and concatenate them together to create a single matrix. The data retains its original statMatrix organization in the variable called 'ensembleMatrix' and is paired with a column ID vector containing the corresponding unit names 'ensembleMatrixColIDs' for logical indexing to select or remove units.

************************************************************************
# List of Required Functions/Toolboxes
************************************************************************
---Plexon Offline Files SDK---
Toolbox created by Plexon to analyze .plx files in Matlab.

---SummarizePLXabbr_BOS.m---
Used when extracting behavioral data from original .plx session files, specifically tailored for use with the .plx files recorded by NJF in Boston. It will not work with files recorded at UCI as the event names used to identify trial boundaries are different.

---SummarizePLXabbr.m---
Used to extract behavioral data from original .plx session files, specifically tailored for use with .plx files recorded at UCI. It will not work on files recorded at Boston as the event names used to identify trial boundaries are different.

---PhaseFreqDetectAbbr.m---
Used for bandpass filtering and hilbert transformation of the LFP data.
