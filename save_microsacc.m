close all

%number of trials
numOfTrials=48; %should be 48 in normal cases

%Read trial details from file
filename = './results/mofo2.txt';
delimiterIn = '\t';
%number of lines to read
headerlinesIn = numOfTrials+1;
%get the lines from file
lines = importdata(filename,delimiterIn,headerlinesIn);

%microsaccades across all trials
final_result={};

%total number of microsaccades
totalMicrosacc=0;
%microsaccade rate over all trials
microsaccRate=[];
%mean microsaccade rate
meanMicrosaccRate=0;


%% just for test (main sequences)
%collect number of main sequences
main_sequences=[];
for i=2:numOfTrials+1
    perLine=strsplit(lines{i,1},'\t'); %Each line
    if perLine{2}=='1'
        main_sequences=[main_sequences,i-1];
    end
end


%%

address=["./results/mofo2.edf";"./results/mofo3.edf"];

total_map={};
total_tmpx={};
total_tmpy={};
total_pupil_data={};

for i=1:size(address,1)
    edfFile=[char(address(i,1))];
    myFile=Edf2Mat(edfFile);
    
    events=myFile.Events.Messages;   %find events of eyelink
    
    %contains events.time correspond to TRIALID n
    %According to number of trials sometimes instead of 24 the array size should be 12 or 8(for each sub-trial)
    startTrial=zeros(1,24);
    %contains events.time correspond to trial OK
    %According to number of trials sometimes instead of 24 the array size should be 12 or 8(for each sub-trial)
    endTrial=zeros(1,24);
    
    samples=myFile.Samples.time;  %get the time each sample was recorded
    
    %collecting start time of each trial
    %("TRIALID n" in events.info corresponds to a time in events.time)
    startTrial=events.time(1,find(contains(events.info,'TRIALID')));
    
    
    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in startTrial
    start_indx_tmp=arrayfun(@(x) find(x==samples,1),startTrial,'un',0);
    
    
    %press_indx=cell2mat(press_indx_tmp);
    % convert cell array which is output of ('un',0) to a vector
    startIndx = cell2mat(start_indx_tmp);
    
    
    %collecting end time of each trial
    %("trial OK" in events.info corresponds to a time in events.time)
    endTrial=events.time(1,find(contains(events.info,'trial OK')));
    
    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in endTrial
    end_indx_tmp=arrayfun(@(x) find(x==samples,1),endTrial,'un',0);
    
    
    %convert cell array which is output of ('un',0) to a vector
    endIndx = cell2mat(end_indx_tmp);
    
    
    SAMPLING=500; %As it is recommended(sample rate)
    TYPE=2;  %velocity type: TYPE=2 recommended
    
    span=100; %It should be changed to 200 for sampling=1000 and 100 for sampling=500
    
    for j=1:length(endTrial)
        %check if the trial is main trial
        if (any(main_sequences(:)==j)&(i==1)|any(main_sequences(:)==j+24)&(i==2))
            %accepted microsaccades of this trial
            total_result=[];
            %trial time record
            trial_time_rec=myFile.Samples.time(startIndx(j):endIndx(j));
            %mapping of the time of each approved data (before elimination, after elimination)
            event_time_map=[1:length(trial_time_rec)];
            %get posx of both left and right eye during each trial
            tmpx=myFile.Samples.posX(startIndx(j):endIndx(j),:);
            
            %get posy of both left and right eye during each trial
            tmpy=myFile.Samples.posY(startIndx(j):endIndx(j),:);
            
            %get the pupil size of the samples in the interval
            pupilSize = myFile.Samples.pupilSize(startIndx(j):endIndx(j),:);
            
            %pupilSize arr counter(remove samples which the increase and decrease of
            %the pupil size is more than 20 units and 200 samples before and after that)
            k=1;
            
            areaThreshold=20; %should be 10 for Sampling=1000 and should be 20 for Sampling=500
            while k<=(size(pupilSize,1)-1)
                if(abs(pupilSize(k+1,1)-pupilSize(k,1))>areaThreshold||abs(pupilSize(k+1,2)-pupilSize(k,2))>areaThreshold)
                    %window off size 400 is flexible according to length of the array
                    if k-span>0 && k+span<size(tmpx,1)
                        pupilSize(k-span:k+span,:) = -1;
                        %the last moment before elimination
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        %the first moment after elimination
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    elseif k+span>size(tmpx,1)
                        pupilSize(k-span:size(pupilSize,1),:) = -1;
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        after_eliminate=trial_time_rec(length(trial_time_rec))-trial_time_rec(1);
                        k=size(tmpx,1)+1;
                    elseif k-span<0
                        %fprintf("K and size(tmpx): %d\t%d\n",k,size(tmpx));
                        pupilSize(1:k+span,:) = -1;
                        before_eliminate=0;
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    end
                    %fprintf("Pupilsize Trial Number %d/ round %d/ Before Eliminate %d/ After Eliminate %d\n",j,k, before_eliminate,after_eliminate);
                    event_time_map(before_eliminate:after_eliminate)= -1;
                else
                    k=k+1;
                end
            end
            
            %tmpx,tmpy counter
            k=1;
            %remove samples in which the data is not recorded due to the
            %subject blink and 200 samples before and after that
            while k<=size(tmpx,1)
                % the window of size 400 is flexible according to length of the
                % array
                if(isnan(tmpx(k,1))||isnan(tmpx(k,2)))
                    if k-span>0 && k+span<size(tmpx,1)
                        tmpx(k-span:k+span,:) = -1;
                        tmpy(k-span:k+span,:) = -1;
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    elseif k+span>size(tmpx,1) && k-span>0
                        tmpx(k-span:size(tmpx,1),:) = -1;
                        tmpy(k-span:size(tmpy,1),:) = -1;
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        after_eliminate=trial_time_rec(length(trial_time_rec))-trial_time_rec(1);
                        k=size(tmpx,1)+1;
                    elseif k+span>size(tmpx,1) && k-span<0
                        tmpx(1:size(tmpx,1),:) = -1;
                        tmpy(1:size(tmpy,1),:) = -1;
                        before_eliminate=0;
                        after_eliminate=trial_time_rec(length(trial_time_rec))-trial_time_rec(1);
                        k=size(tmpx,1)+1;
                        break;
                    elseif k-span<0
                        tmpx(1:k+span,:) = -1;
                        tmpy(1:k+span,:) = -1;
                        before_eliminate=0;
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    end
                    
                    %fprintf("Blink Trial Number %d/ round %d/ Before Eliminate %d/ After Eliminate %d\n",j,k, before_eliminate,after_eliminate);
                    event_time_map(before_eliminate:after_eliminate)= -1;
                else
                    k=k+1;
                end
                %counter=counter+1;
            end
            
            k=1;
            while k<size(event_time_map,2)
                if (event_time_map(k)==-1)
                    %fprintf("time: %d\n",k);
                    event_time_map(k)=[];
                else
                    k=k+1;
                end
            end
            k=1;
            %removing rows whith values -1 in tmpx,tmpy, and pupilSize
            while k<=size(tmpx,1)
                if (tmpx(k,:)==-1 | tmpy(k,:)==-1 | pupilSize(k,:)==-1)
                    tmpx(k,:)=[];
                    tmpy(k,:)=[];
                    pupilSize(k,:)=[];
                else
                    k=k+1;
                end
            end
            %fprintf("length tmpx %d/ length tmpy %d/ length pupilSize %d/ length event_tmp %d\n",length(tmpx),length(tmpy),length(pupilSize),length(event_time_map));
            %This variable is used to calculate miccrosaccade rate (rate=Num/window size(s))
            trialMicrosacc=0;
            
            if(~isempty(tmpx) && ~isempty(tmpy))
                %Engbert algorithm
                if(size(tmpx,1)>=6&&size(tmpy,1)>=6)
                    xx_left=[tmpx(:,1),tmpy(:,1)];
                    v1 = vecvel(xx_left,SAMPLING,TYPE);
                    
                    xx_right=[tmpx(:,2),tmpy(:,2)];
                    v2 = vecvel(xx_right,SAMPLING,TYPE);
                    
                    %               vfac         - [double] velocity factor ("lambda") to determine
                    %                  the velocity threshold for saccade detection
                    %                  (cf. Engbert & Mergenthaler, 2006)
                    %   mindur       - [integer] minimum saccade duration (in samples)
                    %                  (cf. Engbert & Mergenthaler, 2006)
                    
                    VFAC=3;%it was 6 in the main paper [Engbert];
                    MINDUR=3; %it was 3 in the main paper [Engbert];
                    [sac_left, radius_left] = microsacc(xx_left,v1,VFAC,MINDUR);
                    [sac_right, radius_right]=microsacc(xx_right,v2,VFAC,MINDUR);
                    %sac is binocular microsaccades
                    [sac, monol, monor] = binsacc(sac_left,sac_right);
                    
                    result=saccpar(sac);
                    
                    %overshoot correction (remove saccades which the time interval between the offset
                    %of first one and onset of second one is less than 20ms)
                    l=2; %variable for traversing the result array
                    while l<=size(result,1)
                        if((result(l,1))-result((l-1),2)<20)
                            result(l,:)=[];
                        else
                            l=l+1;
                        end
                    end
                    
                    %removing microsaccades with magnitude greater than
                    %2degress
                    l=1; %variable for traversing the result array (one degree is 36 pixels in our config)
                    while l<=size(result,1)
                        if (abs(result(l,8))/36) >2 %maximum threshold for microssacade mag is 2degrees
                            disp(l);
                            result(l,:)=[];
                        else
                            l=l+1;
                        end
                        
                    end
                    
                    %adding to the number of microsaccades
                    trialMicrosacc=trialMicrosacc+size(result,1);
                    total_result=[total_result;result];
                end
            end
            final_result{length(final_result)+1}=total_result;
            total_map{length(total_map)+1}=event_time_map;
            total_tmpx{length(total_tmpx)+1}=tmpx;
            total_tmpy{length(total_tmpy)+1}=tmpy;
        end
    end %this if clause is to check main sequences only
end



% save('./matlab_files/mofo_total_map.mat', 'total_map');
% save('./matlab_files/mofo_total_tmpx.mat', 'total_tmpx');
% save('./matlab_files/mofo_total_tmpy.mat', 'total_tmpy');
save('./matlab_files/mofo_final_microsacc.mat', 'final_result');