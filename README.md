# Sequence and pattern detection
 Code to identify sequences and patterns in inter-spike intervals
This code was generated to identify sequences of (near-) identical inter-spike intervals in neuronal firing patterns.  Additional code for a more general pattern analysis program is also in this repository.

Data are assumed to be in the format shown in the ExampleData.bis file.  THe data need to be loaded into the Matlab workspace for further processing.  Assuming that the workspace name of the data is 'ISI', the next step in the analysis can be ...

S = sequence_arrays_v21(ISI,0.01,100,20);

This generates data about the existence of sequences in the data stream.  

Other parameters describing the ISI data can be obtained with the pattern_analysis routine.  The associated parameter_file_new contains example initial parameters.

R = pattern_analysis(ISI,'parameter_file_new.m');

Summary data can then be extracted with ...

E = extract_basic_parameters_v14(S,R,2);

The code contains hopefully meaningful explanatory wording in the initial Help section, and comments along the way.
