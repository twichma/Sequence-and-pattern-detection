function SRC = sequence_arrays_v21(ISI,thr,num_rep,original_total_rge)
% This function searches for recurring sequences within the input data
% stream ISI (interspike intervals in ms). The analysis uses a relative
% criterion, with which every member of a prospective sequence is matched
% against the corresponding original sequence of ISIs if the duration of
% that ISI is at most thr times different (for example, if thr is set to
% 0.01, and the original sequence is [I01 I02] = [20 50], then a
% repetition of the sequence [I11 I12] will be detected if 0.99*20 <= I11
% <= 1.01*20 AND if 0.99*50 <= I12 <= 1.01*50).  To allow proper memory
% allocation, the user needs to specify the maximal expected sequence length
% (original_total_rge). The analysis is first done on the original data
% (in the sequence ISI), and then repeated on num_rep shuffled 
% representations of the ISI data.  The code also contains a helper 
% function to remove duplicates and overlaps between identified sequences.
%
% Example call: 
% >> S = sequence_arrays_v21(ISI,0.01,100,20);
% 
% This call will analyze the input data stream ISI, using a threshold of
% 0.01. Sequences up to 20 ISIs in length will be analyzed.  For
% statistical comparisons, 100 shuffled representations of the data will be
% analyzed as well.
%
% The paper entitled "Basal Ganglia Neurons in Healthy and Parkinsonian 
% Primates Produce Recurring Sequences of Spikes" (A. Galvan & T. Wichmann)
% is based on version 20 of this routine. Version 21 was prepared for 
% GitHub release, containing additional comments for prospective users.

% Written by TW, 2019-2021. 

global C;

% Initialize output variable
SRC = struct();
for i = 1:original_total_rge
    SRC(i).original_data.input.ISI = [];
    SRC(i).original_data.global.S = {[]};
    SRC(i).original_data.in_range.S = {[]};
    for j = 1:num_rep
        SRC(i).stat(j).input.ISI = [];
        SRC(i).stat(j).global.S = {[]};
        SRC(i).stat(j).in_range.S = {[]};
    end
end

SRC(1).original_data.input.ISI = ISI;
SRC(1).original_data.input.thr = thr;
SRC(1).original_data.input.original_total_range = original_total_rge;
SRC(1).original_data.input.num_rep = num_rep;

%% original data
A = ISI.* ones(length(ISI));% array with length(ISI) column. The complete ISI vector is reproduced in every column of A
C = true(length(ISI),1);    % boolean aray of ISI length, marking those ISIs that are still available to be included in sequences (avoids double-counting of ISIs
B = false(size(A));
total_rge = original_total_rge;

A = A ./ A';                % this matrix contains the results of divisions of each member of A by each member of A (A(1,:) contains all possible divisions by ISI(1), A(2,:) contains all possible divisions by ISI(2) and so on). 
B(abs(A-1)<thr)=true;       % this identifies all division results that are below the threshold as 'true' - the problem is now to identify series of true values
N = false(length(B)-total_rge+1,length(B)-total_rge+1,total_rge);
N(:,:,1) = B(1:length(B)-total_rge+1,1:1:length(B)-total_rge+1);    % the N(:,:,1) array contains a copy of the B array
for i = 2:total_rge         % now we go through the entire range of possible sequence lengths
    N(:,:,i) = N(:,:,i-1) & B(i:length(B)-total_rge+i,i:length(B)-total_rge+i); % ... and see whether there are neighboring truth values - this will end up being true if true values are next to one another, and false if they are not
    if isempty(find(sum(N(:,:,i),2) > 1,1)) % if no sequence is detected, the search is called off
        total_rge = i-1;
        break;
    end
end
for i = total_rge:-1:1      % now we start with the longest detected sequences 
    [SRC(i).original_data.global.S,SRC(i).original_data.in_range.S] = remove_duplicates_and_overlaps(i,N(:,:,i));  % ... and check for duplicates and overlaps
end

SRC(1).original_data.global_stat.actual_total_range = total_rge;

%% statistical analysis (shuffled data)
for nr = 1:num_rep
    disp(['Shuffled sequence # ' num2str(nr)]);
    
    % Generate shuffled version of ISI
    ISIr = ISI(randperm(length(ISI)));
    SRC(1).stat(nr).input.ISI = ISIr;
    
    % (Re-) initialize variables
    A = ISIr.* ones(length(ISIr));
    C = true(length(ISIr),1);    % boolean aray of ISI length, marking those ISIs that are still available to be included in sequences (avoids double-counting of ISIs
    B = false(size(A));
    total_rge = original_total_rge;
    
    A = A ./ A';                % this matrix contains the results of divisions of each member of A by each member of A (B(1,:) contains all possible divisions by ISI(1), B(2,:) contains all possible divisions by ISI(2) and so on). 
    B(abs(A-1)<thr)=true;
    N = false(length(B)-total_rge+1,length(B)-total_rge+1,total_rge);
    N(:,:,1) = B(1:length(B)-total_rge+1,1:1:length(B)-total_rge+1);
    for i = 2:total_rge
        N(:,:,i) = N(:,:,i-1) & B(i:length(B)-total_rge+i,i:length(B)-total_rge+i);
        if isempty(find(sum(N(:,:,i),2) > 1,1))
            total_rge = i-1;
            break;
        end
    end
    for i = total_rge:-1:1
        [SRC(i).stat(nr).global.S,SRC(i).stat(nr).in_range.S] = remove_duplicates_and_overlaps(i,N(:,:,i));
    end
    
    SRC(1).stat(nr).global_stat.actual_total_range = total_rge;
end

end

function [H,K] = remove_duplicates_and_overlaps(rge,N1)
% This function goes through the list of provisional hits, checks for
% duplicates, and assigns the final sequence results. rge is the current
% length of sequences being checked, and N1 is the corresponding N array
% which contains the starting indices of sequences that have a length of
% rge.

global C;                                   % C is a global vector with a list of ISIs.  This start off as a string of zeros, and is gradually replaced by 1s as hits are detected  

% initialize variables
D = ones(size(C));                          % this is used for duplication checking within the same length of sequences
num_unique_seq_global = 0;   
num_unique_seq_in_rge = 0;                   
H = {[]};
K = {[]};
hits = find(sum(N1,2) > 1)';
for i = hits                                % go through N hit-by-hit
    h = find(N1(i,i:end))+i-1;              % find indices of (valid) 1s in each hit row. Numbers less than i are irrelevant as they would have already been found during earlier searches of the same range.  
    
    % *** STEP1: check for overlapping sequences (= sequences with starting points that are within another sequence of the same length)
    m = true(length(h),1);
    for j = 1:length(h)-1
        if m(j)
            for l = j+1:length(h)
                if h(l) < h(j)+rge
                    m(l) = false;
                else
                    break;
                end
            end
        end
    end
    h = h(m == true);
        
    k = h;                                  % this serves to process in-range sequence reporting
    
    % *** STEP 2a: check sequence hits for duplicates GLOBALLY (using h and C arrays)
    m = true(length(h),1);
    for j = 1:length(h)                     % checking for all entries in h ...
        if sum(C(h(j):h(j)+rge-1)) < rge    % are all ISIs within the identified 'sequence' still available?
            m(j) = false;                       % if not, set h(j) provisionally to 0
        end
    end
    h = h(m == true);                            % exclude all zero entries in h
    
    % *** STEP 2b: assign output for global array
    if length(h)>1                          % if h is still longer than 1 entry (length statement works faster than isvector function!)
        num_unique_seq_global = num_unique_seq_global+1;  % increment the unique sequence counter by 1
        H{num_unique_seq_global} = h;                     % assign the index of the first ISI of each repetition of this sequence to H
        for j = 1:length(h)                 % for each entry in the h vector ...
            C(h(j):h(j)+rge-1) = false;     % label the ISIs in C as 'taken'
        end
    end
    
    % *** STEP 3a: check sequence hits for duplicates WITHIN RANGE (using k and D arrays)
    m = true(length(k),1);
    for j = 1:length(k)                     % checking for all entries in k ...
        if sum(D(k(j):k(j)+rge-1)) < rge    % are all ISIs within the identified 'sequence' still available?
            m(j) = false;                       % if not, set k(j) provisionally to 0
        end
    end
    k = k(m == true);                            % exclude all zero entries in k
    
    % *** STEP 3b: assign output for in-range array
    if length(k)>1                          % if k is still longer than 1 entry (length statement works faster than isvector function!)
        num_unique_seq_in_rge = num_unique_seq_in_rge+1;  % increment the unique sequence in range counter by 1
        K{num_unique_seq_in_rge} = k;                     % assign the index of the first ISI of each repetition of this sequence to H
        for j = 1:length(k)                 % for each entry in the h vector ...
            D(k(j):k(j)+rge-1) = false;     % label the ISIs in D as 'taken'
        end
    end
end
end