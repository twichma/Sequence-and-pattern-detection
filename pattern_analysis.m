function R = pattern_analysis(ISI,pfname)
% This function is a collection of routines used for basic pattern analysis purposes
% It uses the input ISI stream ISI to calculate summary statistics
% (duration, mean, median, SD, CV etc) of the input data, and analyzes
% autocorrelograms, power spectra, bursts and declerations, as well as
% pauses in discharge. A variety of input parameters are being used,
% provided in the 'pfname' initialization file.  All output is contained in
% the output structure R.
%
% Example call:
% >> R = pattern_analysis(ISI,'parameter_file_new.m');
%
% This call will carry out the complete series of pattern analysis steps
% using the parameters specified in the parameter file
% 'parameter_file_new.m'.
%
% The paper entitled "Basal Ganglia Neurons in Healthy and Parkinsonian 
% Primates Produce Recurring Sequences of Spikes" (A. Galvan & T. Wichmann)
% is based on version '12/18/2020' of this routine. Version '4/4/2022' does
% only include additional comments (no changes to the code itself).
%
% Written by TW, 4/1/2022


[~, fname, ~] = fileparts(pfname);
eval(fname);
R.pattern_analysis_version = '04/01/2022';
R.ISI = ISI;
R.duration = sum(ISI)/1000;                                       % duration in seconds
R.ISI_percentiles = prctile(ISI,[1,2,5,10,25,50,75,90,95,98,99]);
R.ISI_mean = mean(ISI);
R.ISI_SD = std(ISI);
R.ISI_CV = R.ISI_SD/R.ISI_mean;
R.ISI_SEM = R.ISI_SD/sqrt(length(ISI));
R.firing_rate = length(ISI)/R.duration;
R.num_ISIs = length(ISI);
R.autocorrelogram = ISI_AC(ISI,len,num_rep);
R.powerspectrum = PS_H(ISI,res,iPS_ref,iPS_comp,iPS_method);
R.burst_decel = legendy_new5(ISI,min_length_of_burst,min_length_of_deceleration,surprise_cutoff_burst,surprise_cutoff_deceleration,pause_before_burst);
R.pause_detection = find_pauses(ISI,p_dur);

end % of function

%% Autocorrelation
function autocorrelogram = ISI_AC(ISI,len,num_rep)
% This function calculates a true AC, based on the input data stream (no
% conversion to frequencies).  
% ISI = ISI stream in ms
% len = longest lag in autocorrelation (in ms)
% res = resolution of autocorrelogram (in ms)
% sm_factor = moving average parameter (for instance, with res = 0.5 and an

% Initialize autocorrelogram parameters
autocorrelogram.command.len = NaN;
autocorrelogram.command.num_rep = NaN;
autocorrelogram.correlogram = NaN;
autocorrelogram.reference.correlogram = NaN;

% calculate autocorrelograms
AC = AC_calc(ISI,len);

% generation of reference autocorrelograms and power spectra
AC_r = zeros(num_rep,length(AC));
for n = 1:num_rep
    AC_r(n,:) = AC_calc(ISI(randperm(length(ISI))),len);
end

autocorrelogram.command.len = len;
autocorrelogram.correlogram = AC;
autocorrelogram.reference.correlogram = AC_r;
end % of function

function AC = AC_calc(ISI,len)
% Actual autocorrelogram calculation

ISI_s = zeros(ceil(sum(ISI)),1);
ISI_c = round(cumsum(ISI));
ISI_s(ISI_c) = 1;
max_len = find(ISI_c < sum(ISI)-len,1,'last');
AC = zeros(len,1);
len_red = len-1;

for i = 1:max_len
    AC = AC+ISI_s(ISI_c(i):ISI_c(i)+len_red);
end

end % of function

%% Powerspectrum
function PS = PS_H(ISI,res,iPS_ref,iPS_comp,iPS_method)
% Simple power spectral density calculation - this method replicates that
% by Halliday.

PS.freq = NaN;
PS.powerspectrum = NaN;
PS.conf = NaN;
PS.iPS = [NaN NaN NaN NaN NaN NaN];
D = zeros(ceil(sum(ISI))-1,1);
D(round(cumsum(ISI))) = 1;
D = D-mean(D);
nfft = 2^(ceil(log2(1000/res)));
[P,f,Pc] = pwelch(D,hamming(nfft),0,nfft,1000,'ConfidenceLevel',0.95);
PS.freq = f;
PS.powerspectrum = P;
PS.conf = Pc;
if strcmp(iPS_method,'logarithmic')        % integrates logarithmically
    P = log10(P);
end
PS_ref = sum(P(f(:,1) >= iPS_ref(1) & f(:,1) < iPS_ref(2)));
for n = 1:size(iPS_comp,1)
    PS.iPS(n) = sum(P(f(:,1) >= iPS_comp(n,1) & f(:,1) < iPS_comp(n,2)))/PS_ref;
end
end % of function

%% Burst/deceleration analysis
function R = legendy_new5(ISI,min_length_of_burst,min_length_of_deceleration,surprise_cutoff_burst,surprise_cutoff_deceleration,pause_before_burst)
% This function carries out a modified version of the burst and
% deceleration detection algorithm by Legendy and Salcman 
% (J. Neurophysiology 53:926) on the input ISI data stream.  The sampling
% rate is in Hz.  The user can modify the minimal length of a prospective
% burst (min_length_of_burst), the minimal length of decelerations, and the
% surprise cutoff values for either.  In the original paper, this value was
% set to 10 to detect significant spikes.  In this version of the
% algorithm, the supposed 'burst' or 'deceleration' is made in comparison
% to the entire data stream. The pause_before_burst parameter is the factor
% by which ISIs immediately preceding a burst would have to lengthen in
% order to define a rebound burst.  Because rebound bursts are already
% included in the regular burst statistic, the rebound_burst entry into the
% 'burst' structure gives only the indices of the bursts (out of the burst
% listing ("... the nth burst ...") at which rebound bursts were detected.
% A zero entry denotes that there weren't any bursts that fulfilled the
% pause_before_burst criterion.
%
% The function produces a structure R which contains fields describing 
% bursts and decelerations, including the onset and lengths of these events
% (described as indices of the input ISI stream) the highest rate within 
% bursts or decelerations, the average discharge rate inside the event, the
% pre-event baseline rate, as well as the surprise values for individual
% events.  In field 1 of the structure, summary parameters are defined.
%
% Written 9/21-23/2001, 7/24/2003, 11/27/2003, 7/31/2005, 1/4/2007, 
% 7/23/2008, 10/21-22/2008, 10/23/2020 by Thomas Wichmann.

R = struct();
R.burst = [];
R.deceleration = [];
ISI_burst = false(length(ISI),1);
ISI_decel = false(length(ISI),1);

% Initialize burst and deceleration parameters
% Write input parameters to 'command' field
R.command.min_length_of_burst = min_length_of_burst;
R.command.min_length_of_deceleration = min_length_of_deceleration;
R.command.surprise_cutoff_burst = surprise_cutoff_burst;
R.command.surprise_cutoff_deceleration = surprise_cutoff_deceleration;

% burst parameters (arrays - one entry per burst)
R.burst.begin = NaN;
R.burst.end = NaN;
R.burst.burst_ISIs = {[]};
R.burst.surprise = NaN;
R.burst.rate = NaN;
R.burst.duration = NaN;
R.burst.min_ISI = NaN;
R.burst.max_rate = NaN;
R.burst.max_fold_rate = NaN;
R.burst.number_spikes = NaN;
R.burst.rebound_bursts = NaN;   % boolean - true if burst is a rebound burst, i.e., if the ISI prior to the burst is greater than pause_before_burst * baseline rate
R.burst.LTS_bursts = NaN;       % boolean - true if burst is an LTS_burst
R.burst.Rebound_LTS_bursts = NaN;% boolean - true if LTS burst is a rebound burst

% overall burst statistics
R.burst.number_bursts = NaN;
R.burst.mean_spikes_per_burst = NaN;
R.burst.median_spikes_per_burst = NaN;
R.burst.total_spikes_in_bursts = NaN;
R.burst.baseline_rate = NaN;
R.burst.mean_intra_burst_frequency = NaN;
R.burst.median_intra_burst_frequency = NaN;
R.burst.proportion_time_in_bursts = NaN;
R.burst.proportion_spikes_in_bursts = NaN;
R.burst.proportion_ISIs_in_bursts = NaN;
R.burst.number_rebound_bursts = NaN;                % total number of bursts that fulfill rebound conditions
R.burst.proportion_rebound_bursts = NaN;            % proportion of bursts that fulfill rebound conditions
R.burst.number_LTS_bursts = NaN;                    % total number of LTS bursts (with or without rebound)
R.burst.proportion_LTS_bursts = NaN;                % proportion of LTS bursts (with or without rebound feature) among all bursts
R.burst.number_rebound_LTS_bursts = NaN;            % total number of LTS bursts that fulfill rebound fieature
R.burst.proportion_rebound_LTS_bursts = NaN;        % proportion of 'rebound LTS' bursts among all bursts
R.burst.proportion_LTS_rebound_bursts = NaN;        % proportion of 'rebourn LTS' bursts among all rebound bursts

% deceleration parameters (arrays - one entry per burst)
R.deceleration.begin = NaN;
R.deceleration.end = NaN;
R.deceleration.deceleration_ISIs = {[]};
R.deceleration.surprise = NaN;
R.deceleration.rate = NaN;
R.deceleration.duration = NaN;
R.deceleration.max_ISI = NaN;
R.deceleration.min_rate = NaN;
R.deceleration.min_fold_rate = NaN;
R.deceleration.number_spikes = NaN;

% overall deceleration statistics
R.deceleration.number_decelerations = NaN;
R.deceleration.mean_spikes_per_deceleration = NaN;
R.deceleration.median_spikes_per_deceleration = NaN;
R.deceleration.total_spikes_in_decelerations = NaN;
R.deceleration.baseline_rate = NaN;
R.deceleration.mean_intra_deceleration_frequency = NaN;
R.deceleration.median_intra_deceleration_frequency = NaN;
R.deceleration.proportion_time_in_decelerations = NaN;
R.deceleration.proportion_spikes_in_decelerations = NaN;

% unmodulated parameters
R.unmodulated.ISI_boolean = NaN;
R.unmodulated.ISI = NaN;
R.unmodulated.duration = NaN;
R.unmodulated.ISI_percentiles = NaN;
R.unmodulated.ISI_mean = NaN;
R.unmodulated.ISI_SD = NaN;
R.unmodulated.ISI_CV = NaN;
R.unmodulated.ISI_SEM = NaN;
R.unmodulated.firing_rate = NaN;
R.unmodulated.number_ISIs = NaN;
R.unmodulated.proportion_unmodulated_spikes = NaN;
R.unmodulated.proportion_unmodulated_time = NaN;

% Analysis of bursts
s = brute_force_surprise_v2(ISI,1,R.command.surprise_cutoff_burst,R.command.min_length_of_burst,25);
R.burst.begin = s.start;
R.burst.end = s.end;
R.burst.number_spikes = s.len_sp;
R.burst.mean_spikes_per_burst = mean(s.len_sp);
R.burst.median_spikes_per_burst = median(s.len_sp);
R.burst.total_spikes_in_bursts = sum(s.len_sp);
R.burst.mean_intra_burst_frequency = mean(s.rate);
R.burst.median_intra_burst_frequency = median(s.rate);
R.burst.duration = s.dur';
R.burst.surprise = s.surp;
R.burst.rate = s.rate;
R.burst.max_rate = 1000./s.min';
R.burst.min_ISI = s.min';
R.burst.baseline_rate = sum(ISI)/length(ISI);
R.burst.max_fold_rate = s.max_fold_rate';
R.burst.number_bursts = length(R.burst.begin(~isnan(R.burst.begin)));
R.burst.proportion_spikes_in_bursts = s.prop_sp;
R.burst.proportion_ISIs_in_bursts = s.prop_ISI;
R.burst.proportion_time_in_bursts = s.prop_t;

% 'Rebound' burst analysis
rebound_bursts = false(length(s.start),1);
for burst_num = 1:length(s.start)
    if s.start(burst_num) > 1 && (ISI(s.start(burst_num)-1) > pause_before_burst*R.burst.baseline_rate)
        rebound_bursts(burst_num) = true;
    end
    R.burst.burst_ISIs{burst_num} = ISI(s.start(burst_num):s.end(burst_num));
    ISI_burst(s.start(burst_num):s.end(burst_num)) = true;
end
R.burst.rebound_bursts = rebound_bursts;
R.burst.number_rebound_bursts = sum(R.burst.rebound_bursts);
R.burst.proportion_rebound_bursts = R.burst.number_rebound_bursts/R.burst.number_bursts;

% LTS burst analysis
R.burst.number_bursts = length(R.burst.begin(~isnan(R.burst.begin)));
R.burst.LTS_bursts = false(R.burst.number_bursts,1);
R.burst.Rebound_LTS_bursts = false(R.burst.number_bursts,1);

for i = 1:R.burst.number_bursts             % for all available bursts ...
    if diff(R.burst.burst_ISIs{i}) > 0                        % this checks whether ALL diff(ISIs) are > 0, i.e., whether the ISIs are uniformly increasing!  
        R.burst.LTS_bursts(i) = true;
        if R.burst.rebound_bursts(i)                  % this assigns the reboud_LTS burst label in case an LTS burst is also a reoubd burst  
            R.burst.Rebound_LTS_bursts(i) = true;
        end
    end
end

% assign entries in the output structure
R.burst.number_LTS_bursts = sum(R.burst.LTS_bursts);               % number of LTS bursts detected
R.burst.proportion_LTS_bursts = R.burst.number_LTS_bursts/R.burst.number_bursts;     % proportion of bursts that fulfill the monotonous increase criterion
R.burst.number_rebound_LTS_bursts = sum(R.burst.Rebound_LTS_bursts);
R.burst.proportion_rebound_LTS_bursts = R.burst.number_rebound_LTS_bursts/R.burst.number_bursts;  % proportion of bursts that fulfill rebound & LTS criteria
R.burst.proportion_LTS_rebound_bursts = R.burst.number_rebound_LTS_bursts/sum(R.burst.rebound_bursts);   % proportion of rebound bursts that fulfill LTS criteria

% Analysis of decelerations
s = brute_force_surprise_v2(ISI,0,R.command.surprise_cutoff_deceleration,R.command.min_length_of_deceleration,25);
for dec_num = 1:length(s.start)
    R.deceleration.deceleration_ISIs{dec_num} = ISI(s.start(dec_num):s.end(dec_num));
    ISI_decel(s.start(dec_num):s.end(dec_num)) = true;
end
R.deceleration.begin = s.start;
R.deceleration.end = s.end;
R.deceleration.number_spikes = s.len_sp;
R.deceleration.mean_spikes_per_deceleration = mean(s.len_sp);
R.deceleration.median_spikes_per_deceleration = median(s.len_sp);
R.deceleration.total_spikes_in_decelerations = sum(s.len_sp);
R.deceleration.mean_intra_deceleration_frequency = mean(s.rate);
R.deceleration.median_intra_deceleration_frequency = median(s.rate);
R.deceleration.proportion_time_in_decelerations = s.prop_t;
R.deceleration.proportion_spikes_in_decelerations = s.prop_sp;
R.deceleration.duration = s.dur';
R.deceleration.surprise = s.surp;
R.deceleration.rate = s.rate;
R.deceleration.min_rate = 1000./s.max';
R.deceleration.max_ISI = s.max';
R.deceleration.baseline_rate = sum(ISI)/length(ISI);
R.deceleration.min_fold_rate = s.min_fold_rate';
R.deceleration.number_decelerations = length(R.deceleration.begin(~isnan(R.deceleration.begin)));

% Analysis of unmodulated data
ISI_unmod = ~(ISI_burst | ISI_decel);
R.unmodulated.ISI_boolean = ISI_unmod;
R.unmodulated.ISI = ISI(ISI_unmod);
R.unmodulated.duration = sum(R.unmodulated.ISI)/1000;                                       % duration in seconds
R.unmodulated.ISI_percentiles = prctile(R.unmodulated.ISI,[1,2,5,10,25,50,75,90,95,98,99]);
R.unmodulated.ISI_mean = mean(R.unmodulated.ISI);
R.unmodulated.ISI_SD = std(R.unmodulated.ISI);
R.unmodulated.ISI_CV = R.unmodulated.ISI_SD/R.unmodulated.ISI_mean;
R.unmodulated.ISI_SEM = R.unmodulated.ISI_SD/sqrt(length(R.unmodulated.ISI));
R.unmodulated.firing_rate = length(R.unmodulated.ISI)/R.unmodulated.duration;
R.unmodulated.number_ISIs = length(R.unmodulated.ISI);
R.unmodulated.proportion_unmodulated_spikes = R.unmodulated.number_ISIs/length(ISI);
R.unmodulated.proportion_unmodulated_time = sum(R.unmodulated.ISI)/sum(ISI);

end % of function

function s = brute_force_surprise_v2(ISI,mode,sco,min_len,max_len)
% A version of the surprise calculation that takes all possible surprise
% values into account.
% mode = 0 is for decelerations, more 1 is for bursts
% sco = surprise cut-off value
% min_len = minimal burst/decelration length
% max_len = maximal eveluation length of individual bursts (can be
% arbitratily set to 25)

% Initialize variables
ISI_comp = mean(ISI);                           % the average ISI is set as the ISI for comparison
ISI_tf = false(length(ISI),1);                  % ISI_tf is a boolean vector of the same length as ISI
s.surp = NaN(length(ISI),1);                    % each ISI can have a surprise value
s.start = s.surp;                               % burst start
s.end = s.surp;                                 % burst end
s.dur = NaN;
s.max = NaN;
s.min = NaN;
s.rate = NaN;
s.max_fold_rate = NaN;
s.min_fold_rate = NaN;
s.prop_ISI = NaN;
s.prop_sp = NaN;
s.prop_t = NaN;
s.len_ISI = NaN;
s.len_sp = NaN;

for i = 1:length(ISI)-max_len                   % we are looking through all ISIs that allow max_len bursts to occur 
    C = cumsum(ISI(i:i+max_len))';              % for each spike, we are generating the cumulative sum, going out to max_length ISIs later
    ss = -log((2*mode-1)*(mode-poisscdf(1:length(C),C(1:end)/ISI_comp))); % surprise series, comparing the number of ISIs found, with the number of ISIs expected
    [surp,ind] = max(ss(min_len:end));                         % surprise series, comparing ISIs found, with ISIs expected
    if surp >= sco                              % if the maximal surprise value was higher than the cut-off
        s.start(i) = i;                         % provisional definition of start... (note: this is set to i!)
        s.end(i) = i+min_len-1+ind-1;           % ... and end of burst
        s.surp(i) = surp;                       % ... and setting the surprise value
    end
end
s.surp = s.surp(~isnan(s.surp));                % we find the actual entried into the s.surp, s.start, and s.end arrays
s.start = s.start(~isnan(s.start));
s.end = s.end(~isnan(s.end));

% Removing overlapping hi-surprise components of the ISI stream
[surp_sorted,ind] = sort(s.surp,'descend');     % and sort them
for i = 1:length(surp_sorted)
    if ~any(ISI_tf(s.start(ind(i)):s.end(ind(i))))  % if all entries into ISI_tf that would correspond to the burst are still zero
        ISI_tf(s.start(ind(i)):s.end(ind(i))) = true;   % the ISI_tf values will be set to true 
    else                                        % if some or all of the ISIs are already taken up ...
        s.surp(ind(i)) = NaN;                   % we blot the corresponding burst out
        s.start(ind(i)) = NaN;
        s.end(ind(i)) = NaN;
    end
end

s.surp = s.surp(~isnan(s.surp));                % 
s.start = s.start(~isnan(s.start));
s.end = s.end(~isnan(s.end));
% note that s.start, s.end, and s.surp are all indices referring to ISIs

% set s_len in two different flavors
s.len_ISI = s.end-s.start+1;
s.len_sp = s.len_ISI + 1;

% calculate duration, highest and lowest ISI values within identified
% surprising sequences

for i = 1:length(s.start)
    s.dur(i) = sum(ISI(s.start(i):s.end(i)));
    s.max(i) = max(ISI(s.start(i):s.end(i)));
    s.min(i) = min(ISI(s.start(i):s.end(i)));
end

% calculate other parameters
if ~isempty(s.start)
    s.rate = 1000*s.len_ISI./s.dur';
    s.max_fold_rate = mean(ISI)./s.min;
    s.min_fold_rate = mean(ISI)./s.max;
    s.prop_ISI = sum(ISI_tf)/length(ISI);
    s.prop_sp = (sum(ISI_tf)+length(s.start))/length(ISI);
    s.prop_t = sum(ISI(ISI_tf))/sum(ISI);
end
end

%% *** pause analysis
function P = find_pauses(A,p_dur)
% as opposed to the legendy routine, this program finds true pauses, i.e.,
% cessations of neuronal activity for a set number of ms (p_dur)

% Initialize pause detection parameters
P.pause.command.time = clock;
P.pause.command.p_dur = p_dur;
P.pause.num_spikes = NaN;
P.pause.begin = NaN;
P.pause.ISIs{1} = NaN;
P.pause.duration = NaN;
P.pause.mean_ISI_in_pauses = NaN;
P.pause.SD_ISI_in_pauses = NaN;
P.pause.SEM_ISI_in_pauses = NaN;
P.pause.percentiles_ISI_in_pauses = NaN;
P.pause.proportion_spikes_in_pauses = NaN;
P.pause.proportion_time_in_pauses = NaN;
P.pause.mean_ISI_in_single_pause = NaN;

% Pause detection
short_ISI = find(A < p_dur);
long_ISI = find(A(1:short_ISI(end)) >= p_dur);                        % identify long ISIs. In order to define a pause, it needs to be followed by at least one short ISI.
pause_start = long_ISI(diff([-1;long_ISI]) > 1);    % identify those that are separate by at least one ISI
if ~isempty(pause_start)
    pause_length = zeros(size(pause_start));            % initialize array to identify the length of the pause
    P.pause.num_spikes = zeros(size(pause_start));
    P.pause.begin = pause_start';
    P.pause.duration = zeros(size(pause_start));
    P.pause.mean_ISI_in_single_pause = zeros(size(pause_start));
    P.pause.ISIs = {};
    for i = 1:length(pause_start)                       % for all pauses
        pause_length(i) = short_ISI(find(short_ISI > pause_start(i),1,'first'))-pause_start(i);
        P.pause.num_spikes(i) = pause_length(i);
        P.pause.ISIs{i} = A(pause_start(i):pause_start(i)+pause_length(i));
        P.pause.duration(i) = sum(P.pause.ISIs{i});
        P.pause.mean_ISI_in_single_pause(i) = mean(P.pause.ISIs{i});
    end
    P.pause.mean_ISI_in_pauses = mean(P.pause.mean_ISI_in_single_pause);
    P.pause.SD_ISI_in_pauses = std(P.pause.mean_ISI_in_single_pause);
    P.pause.SEM_ISI_in_pauses = P.pause.SD_ISI_in_pauses/sqrt(length(pause_start));
    P.pause.proportion_spikes_in_pauses = sum(pause_length)/length(A);
    P.pause.proportion_time_in_pauses = sum(P.pause.duration)/sum(A);
end
end % of function