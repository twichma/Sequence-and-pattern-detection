% **** burst/deceleration detection ****
min_length_of_burst = 2;            % burst detection - minimal number of ISIs in burst
min_length_of_deceleration = 1;     % deceleration detection - minimal number of ISIs in deceleration
surprise_cutoff_burst = 3;          % burst detection - surprise cutoff value
surprise_cutoff_deceleration = 3;   % deceleration detection - surprise cutoff value
pause_before_burst = 2;             % denotes how slow a cell needs to get right before a burst, so that the burst is classified as a rebound burst (e.g., '3' would mean that the ISI before the burst would have to be three times the mean ISI of the cell.

% **** pause detection ****
p_dur = 500;                        % pause detection - specifies the minimal length (in ms) of the pause; if p_dur is defined (not []), p_rel is not considered.

% **** Autocorrelation, powerspectrum routine ****
len  = 1000;                        % length of autocorrelation lag
num_rep = 250;                      % used for statistics - number of shuffled data representations

% **** Power spectrum routine ****
res = 1;                            % determines the 'resolution' - time base of the input ISIs.  In the usual case of an input resolution of 1 ms/sample, res is 1.
iPS_ref = [1 100];                  % this are the reference boundaries for the integrated power spectrum (e.g., iPS_ref = [1 100]; would refer to reference boundaries of 1 Hz (lower limit) and 100 Hz (upper limit))
iPS_comp = [1 3;3 8;8 13;13 30;30 100]; % these are the component areas for the integrated power spectra (e.g., iPS_comp = [1 3;3 8] would result in integration of 1->3 Hz and 3->8 Hz powerspectra
iPS_method = 'linear';               % this string defines whether the linear ('linear') or logarithmic ('logarithmic') power spectra are integrated