function [corr] = test_correlation2(usacc_on,pr,trial_props,trialRows)
%% window backward will go back from time of button press until the start of the trial at the most. 
%% Window forward will go forward into the next trial and till the button press in the next trial at the most. 
u = usacc_on; %list of microsaccade onsets
num_events = size(pr,1); %no. of button presses. 'pr' is a list of indices of button presses
win_back = 2000; %no. of samples to go backward in time from instant of button press 2000 samples = 4000ms    
win_for = 1000; %no. of samples to go forward in time from instant of button press.
sumdat = zeros(win_back + win_for,1);
total_window = sumdat;
trial_starts = trial_props(:,1);%list of start indices of trials
trial_stops = trial_props(:,2);%list of stop indices of trials
corrData = zeros(num_events,5); %for my reference to keep track of the data
%note that for this experiment, a button was pressed only once during a trial and once the button was pressed, the trial ended
for event = 1:num_events %iterate through the button presses
    start = trial_starts(event); %get start index of current trial
    stop = trial_stops(event); %get stop index of current trial
    curr_event = pr(event); %get index of button press
    corrStart = max(start, curr_event - win_back);%start index of data that must be considered for correlation 
    % - either the start of the trial or 4s back into the trial
    if event ~= num_events %if it isn't the last event
        corrStop = min(curr_event + win_for,pr(event+1));%stop index of data that must be considered for correlation 
        % - either 2s forward from the button press or until the next button press (next trial)
    else %if it is the last event
        corrStop = stop;
    end
    corrData(event,1) = curr_event;
    corrData(event,[2 3 5]) = [start stop corrStop];
    startStep = win_back - (curr_event - corrStart);%row at which the data starts in the 3000 step window
    stopStep = win_back + corrStop - curr_event;%row at which the data stops in the 3000 step window
    to_add = zeros(win_back + win_for,1); %stores all the indices during a trial that must be considered for correlation
    corrData(event,7) = stopStep;
    if startStep == 0 %if the correlation starts at the start of the trial
        startStep = 1;
        to_add(startStep:stopStep) = u(corrStart+1:corrStop); %get the relevant portion of the microsaccade onset data
        corrData(event,4) = corrStart+1;
    else
        to_add(startStep:stopStep) = u(corrStart:corrStop); %get the relevant portion of the microsaccade onset data
        corrData(event,4) = corrStart;
    end
    corrData(event,6) = startStep;
    sumdat = sumdat + to_add; %cumulatively add all the relevant portions of each trial
    to_add = zeros(win_back + win_for,1);
    to_add(startStep:stopStep) = 1;
    total_window = total_window + to_add; %keeping a count of the number of windows of data added
end
    
if ~isempty(trialRows) %applicable only if you are picking out certain trials to perform correlation
    corrData = corrData(trialRows,:);
    sumdat = zeros(win_back + win_for,1);
    total_window = sumdat;
    for trial = 1:size(corrData,1)
        to_add = zeros(win_back + win_for,1);
        to_add(corrData(trial,6):corrData(trial,7)) = u(corrData(trial,4):corrData(trial,5));
        sumdat = sumdat + to_add;
        to_add = zeros(win_back + win_for,1);
        to_add(corrData(trial,6):corrData(trial,7)) = 1;
        total_window = total_window + to_add;
    end 
end  
corr = sumdat./total_window * 500; 
%divide the cumulative periods of microsaccade onsets by the total number of windows converted into total time by multiplying by the sample rate (500 Hz)
%This will give us the correlations of microsaccade onsets to button presses in the form of a microsaccade rate