function experiment(sId)

% Here we call some default settings for Psychtoolbox
PsychDefaultSetup(2);
%PsychDebugWindowConfiguration

% Initialization
whichScreen = 0 ; %allow to choose the display if there's more than one
white = WhiteIndex(whichScreen);
[w, rect] = PsychImaging('OpenWindow', whichScreen, white);

[xCenter, yCenter]=RectCenter(rect); %find out the center of the screen
Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); %window transparency


%Loading images
nonIllusoryImg=imread('./images/single_snake_non_illusory.bmp');
illusoryImg=imread('./images/single_snake_illusory.bmp');


%position and color of the centered red dot
dotXpos = xCenter;
dotYpos = yCenter;
dotColor = [1 0 0];

% Dot size in pixels 0.25 degree of visual angle
% in our config one degree is 36 pixels
dotSizePix = 9;

%############################
%% Eyelink: start
% dummymode=0;
% if ~EyelinkInit(dummymode, 1)
%     fprintf('Eyelink Init aborted.\n');
%     Eyelink('Shutdown');
%     return;
% end
% eye_used = -1;
% el=EyelinkInitDefaults(w);
%   % open file to record data to
% edfFile=[sId '.edf'];
% status= Eyelink('openfile',edfFile);
% if status~= 0
%     error('openfile error, status: ', status);
% end
% Eyelink('trackersetup');
% Eyelink('startrecording');
% WaitSecs(1);
% Eyelink('Message', 'SYNC  TIME');
%
% disp('1');
%
% if Eyelink( 'NewFloatSampleAvailable') > 0
%     evt = Eyelink( 'NewestFloatSample');
%     if eye_used ~= -1
%         x = evt.gx(eye_used+1);
%         y = evt.gy(eye_used+1);
%     else % if we don't, first find eye that's   being tracked
%         eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
%         %         if eye_used == el.BINOCULAR; % if both eyes are tracked
%         %             eye_used = el.RIGHT_EYE; % use left eye
%         %         end
%     end
% end
%% Eyelink: end
%############################
%Design 48 trials half of them illusory condition(1) and the other half
%non-illusory condition(0)

%numberOfTrials
numOfTrials=2;
%permutation of trials
r = [ones(1, (numOfTrials)/2 - 1), zeros(1, numOfTrials/2)];
A(2:numOfTrials)=r(randperm(numOfTrials-1));
A(1)=1;  %force illusion to be the first trial
% The avaliable keys to press
escapeKey = KbName('ESCAPE');
%Answer trials
spaceKey = KbName('SPACE');


%create baseRect
width=288;
height=288;
baseRect=[0 0 width height];


%make image into texture
texture1 = Screen('MakeTexture', w, illusoryImg);
texture2 = Screen('MakeTexture', w, nonIllusoryImg);

%open a file with subject's name in order to save relevant information
fileName=[sId,'.txt'];
fid=fopen(fileName,'a+');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','trial_num','condition','directionOfMotion','angle','pressedTimes','releasedTimes','RealTransitionsToStationary','RealTransitionsToRotation');

% the exact time of transition towards rotation state in control condition
exactRotationTime=[];
% the exact time of transition towards stationary state in control condition
exactQuitTime=[];


for i=1:length(A)
    pressCnt=1; %number of pressTime to send message to eyelink
    releaseCnt=1; %number of releaseTime to send message to eyelink
    
    %############################
    %%    Start Trial:start
    %     Eyelink('Message','TRIALID %d', i);
    %     mes=['start trial' num2str(i)];
    %     Eyelink('command','resord_status_messasge "%s"',mes);
    %%    Start Trial:end
    %############################
    
    Screen('DrawDots', w, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
    Screen('Flip',w);
    
    
    % Now we have drawn to the screen we wait for a keyboard button press (any
    % key) to start the trial.
    KbStrokeWait;
    
    %decide for rotation direction (1 clockwise, 0 counter clockwise)
    clockwise=rem(i,2);
    
    
    %draw base rect and replace them around the compass randomly (0,30,60,90,120,150)
    r=randi([1 6]);
    
    dstRects = nan(4, 2);
    
    if (r==1) %0 degree
        yPos = yCenter;
        xPos = [xCenter-324 xCenter+324];
        for j = 1:2
            dstRects(:, j) = CenterRectOnPointd(baseRect, xPos(j), yPos);
        end
        
    elseif(r==2) %90 degrees
        yPos = [yCenter+324 yCenter-324];
        xPos = xCenter;
        for j = 1:2
            dstRects(:, j) = CenterRectOnPointd(baseRect, xPos, yPos(j));
        end
        
    elseif(r==3) %30 degrees
        yPos = [(yCenter-(1/2*324)) (yCenter+(1/2*324))];
        xPos = [(xCenter-(sqrt(3)/2*324)) (xCenter+(sqrt(3)/2*324))];
        for j = 1:2
            dstRects(:, j) = CenterRectOnPointd(baseRect, xPos(j), yPos(j));
        end
        
    elseif(r==4)%60 degrees
        yPos = [(yCenter-(sqrt(3)/2*324)) (yCenter+(sqrt(3)/2*324))];
        xPos = [(xCenter-(1/2*324)) (xCenter+(1/2*324))];
        for j = 1:2
            dstRects(:, j) = CenterRectOnPointd(baseRect, xPos(j), yPos(j));
        end
    elseif(r==5) %120 degrees
        yPos = [(yCenter+(sqrt(3)/2*324)) (yCenter-(sqrt(3)/2*324))];
        xPos = [(xCenter-(1/2*324)) (xCenter+(1/2*324))];
        for j = 1:2
            dstRects(:, j) = CenterRectOnPointd(baseRect, xPos(j), yPos(j));
        end
        
    elseif(r==6) %150 degrees
        yPos = [(yCenter+(1/2*324)) (yCenter-(1/2*324))];
        xPos = [(xCenter-(sqrt(3)/2*324)) (xCenter+(sqrt(3)/2*324))];
        for j = 1:2
            dstRects(:, j) = CenterRectOnPointd(baseRect, xPos(j), yPos(j));
        end
    end
    
    % A matrix used to keep pressTimes and releaseTimes in order to be saved
    % in a file
    KeyTime=[];
    if(A(i)==1) %main condition
        
        %this variable is for indexing pressTimes and holdingTimes
        %arrays
        step=1;
        
        %the arrays below hold the pressTime and holdingTime of the current
        %trial
        pressTimes=zeros(1, 30);
        holdingTimes=zeros(1, 30);
        
        %draw textures and center dot
        Screen('DrawTextures', w, texture1, [],dstRects);
        Screen('DrawDots', w, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
        Screen('Flip', w);
        
        %the trial should last for 30secs
        trialTime=tic;
        
        startSecs = GetSecs;
        while toc(trialTime)<30
            
            [ keyIsDown, timeSecs, keyCode ] = KbCheck;
            keyCode = find(keyCode, 1);
            
            if keyIsDown
                pressTime=0;
                if keyCode == escapeKey
                    break;
                elseif keyCode==spaceKey
                    %############################
                    %%    PressTime:start
                    %     Eyelink('Message','PressTime %d', pressCnt);
                    %     mes=['pressTime' num2str(pressCnt)];
                    %     Eyelink('command','resord_status_messasge "%s"',mes);
                    %%    pressTime:end
                    %############################
                    fprintf('"%s" typed at time %.3f seconds\n', KbName(keyCode), timeSecs - startSecs);
                    pressCnt=pressCnt+1;
                    pressTime= timeSecs - startSecs;
                    if(timeSecs - startSecs>30)
                        break;
                    end
                    
                end
                
                % If the user holds down a key, KbCheck will report multiple events.
                % To condense multiple 'keyDown' events into a single event, we wait until all
                % keys have been released.
                
                %[startTime]=WaitSecs(0);
                startTime=GetSecs;
                %disp(startTime);
                value=30-(startTime-startSecs);
                [secs, keyCode, ~] = KbReleaseWait([],startTime+value); %startTime+value
                
                %############################
                %%    ReleaseTime:start
                %     Eyelink('Message','ReleaseTime %d', releaseCnt);
                %     mes=['ReleaseTime' num2str(releaseCnt)];
                %     Eyelink('command','resord_status_messasge "%s"',mes);
                %%    ReleaseTime:end
                %############################
                
                %save release time as well as hold time
                holdTime= secs - timeSecs;
                releaseTime=pressTime + holdTime;
                fprintf('"%s" released after %.3f seconds\n', "spaceKey",secs - timeSecs);
                releaseCnt=releaseCnt+1;
                pressTimes(1,step)=pressTime;
                holdingTimes(1,step)=holdTime;
                releaseTimes(1,step)=releaseTime;
                step=step+1;
                
                KeyTime=[KeyTime;pressTime,releaseTime];
            end
            
        end
        
        %save trial number and trial condition
        fprintf(fid,'%d\t',i);
        fprintf(fid,'%d\t',A(i));
        
        %save rotation direction in the file
        fprintf(fid,'%d\t',clockwise+1);
        
        %save angle in the file
        fprintf(fid,'%d\t',r);
        
        %saving information of pressedTimes and releasedTimes in text file
        for k = 1:size(KeyTime,1)
            fprintf(fid,'%g ',KeyTime(k,1));
        end
        fprintf(fid,'\t');
        for k = 1:size(KeyTime,1)
            fprintf(fid,'%g ',KeyTime(k,2));
        end
        fprintf(fid,'\n');
        
    elseif(A(i)==0) %control condition
        %used to save time records temporarily
        startRec=[];
        stopRec=[];
        
        % Set the inital rotation angle randomly and in increment per frame to 0.06
        % degrees
        angles = rand(1, 2) .* 360;
        degPerFrame = 0.06;
        
        % Sync us and get a time stamp
        Screen('Flip', w);
        
        
        % Maximum priority level
        topPriorityLevel = MaxPriority(w);
        Priority(topPriorityLevel);
        
        %counter on holdingTime arr
        j=1;
        
        
        % each trial should last about 30 secs
        trialTime=tic;
        %psychtoolbox gets current time
        startSecs = GetSecs;
        %flag to check key stroke
        keyStroke=false;
        pressTime=0;
        while(toc(trialTime)<30)
            %rotation phase
            Screen('DrawDots', w, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
            Screen('DrawTextures', w, texture2, [],dstRects,angles);
            
            % Flip to the screen
            Screen('Flip', w);
            
            % check the press and release patterns from the most recently main
            % condition trial
            if(round(toc(trialTime),1) == round(pressTimes(1,j),1))
                
                % save record to store it as an exact time towards stationary
                % state
                startRec=[startRec,toc(trialTime)];
                %stop phase
                tstart=tic;
                while(toc(tstart)<holdingTimes(1,j))
                    Screen('DrawDots', w, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
                    Screen('DrawTextures', w, texture2, [],dstRects,angles);
                    Screen('Flip', w);
                    
                    [~, timeSecs, keyCode ] = KbCheck;
                    
                    if keyCode(spaceKey)
                        if ~keyStroke
                            keyStroke=true;
                            
                            %############################
                            %%    PressTime:start
                            %     Eyelink('Message','PressTime %d', pressCnt);
                            %     mes=['pressTime' num2str(pressCnt)];
                            %     Eyelink('command','resord_status_messasge "%s"',mes);
                            %%    pressTime:end
                            %############################
                            
                            fprintf('"%s" typed at time %.3f seconds\n', KbName(keyCode), timeSecs - startSecs);
                            pressCnt=pressCnt+1;
                            pressTime= timeSecs - startSecs;
                            
                        end
                    elseif ~keyCode(spaceKey)
                        if keyStroke
                            keyStroke=false;
                            %timeNow=GetSecs;
                            %releaseTime=pressTime+(timeNow-timeSecs);
                            releaseTime=toc(trialTime);
                            %############################
                            %%    ReleaseTime:start
                            %     Eyelink('Message','ReleaseTime %d', releaseCnt);
                            %     mes=['ReleaseTime' num2str(releaseCnt)];
                            %     Eyelink('command','resord_status_messasge "%s"',mes);
                            %%    ReleaseTime:end
                            %############################
                            releaseCnt=releaseCnt+1;
                            fprintf('"%s" released after %.3f seconds\n',"spaceKey" ,releaseTime-pressTime); %save release time as well as hold time
                            KeyTime=[KeyTime;pressTime,releaseTime];
                        end
                    end
                end
                % save record to store it as an exact time towards rotation
                % state
                stopRec=[stopRec,toc(trialTime)];
                j=j+1;
                if(j==31)
                    j=1;
                end
            end
            
            % Increment the angle for the next drawing loop
            if(clockwise==1)
                %clockwise
                angles = angles + degPerFrame;
            elseif(clockwise==0)
                %counter clockwise
                angles = angles - degPerFrame;
            end
            
            [~,secs,keyCode] = KbCheck;
            %the moment when the space bar is released
            %(in case the button is released after the inner while)
            if(~keyCode(spaceKey))
                if keyStroke
                    keyStroke=false;
                    disp("biroon vel kardam :)");
                    %releaseTime=pressTime+(secs-timeSecs);
                    releaseTime=toc(trialTime);
                    %############################
                    %%    ReleaseTime:start
                    %     Eyelink('Message','ReleaseTime %d', releaseCnt);
                    %     mes=['ReleaseTime' num2str(releaseCnt)];
                    %     Eyelink('command','resord_status_messasge "%s"',mes);
                    %%    ReleaseTime:end
                    %############################
                    releaseCnt=releaseCnt+1;
                    fprintf('"%s" released after %.3f seconds\n', "spaceKey",releaseTime-pressTime); %save release time as well as hold time
                    KeyTime=[KeyTime;pressTime,releaseTime];
                end
            elseif keyCode(spaceKey)
                if ~keyStroke
                    keyStroke=true;
                    disp("biroon space zadam :)");
                    %############################
                    %%    PressTime:start
                    %     Eyelink('Message','WrongPress %d', pressCnt);
                    %     mes=['pressTime' num2str(pressCnt)];
                    %     Eyelink('command','resord_status_messasge "%s"',mes);
                    %%    pressTime:end
                    %############################
                    fprintf('"%s" typed at time %.3f seconds\n', "spaceKey", secs - startSecs);
                    pressCnt=pressCnt+1;
                    pressTime= secs - startSecs;
                    
                end
                
            end
            
        end
        
        %save trial number and trial condition
        fprintf(fid,'%d\t',i);
        fprintf(fid,'%d\t',A(i));
        
        %save rotation direction in the file
        fprintf(fid,'%d\t',clockwise+1);
        
        %save angle in the file
        fprintf(fid,'%d\t',r);
        
        %saving information of pressedTimes and releasedTimes in text file
        for k = 1:size(KeyTime,1)
            fprintf(fid,'%g ',KeyTime(k,1));
        end
        fprintf(fid,'\t');
        
        for k = 1:size(KeyTime,1)
            fprintf(fid,'%g ',KeyTime(k,2));
        end
        fprintf(fid,'\t');
        % save current trial transition records
        for k = 1:size(startRec,2)
            fprintf(fid,'%g ',startRec(k));
        end
        fprintf(fid,'\t');
        
        for k = 1:size(stopRec,2)
            fprintf(fid,'%g ',stopRec(k));
        end
        fprintf(fid,'\n');
    end
    
    
    %############################
    %     Eyelink('Message','TRIAL_RESULT', i);
    %     Eyelink('Message','trial OK');
    %############################
    
end
%closing the file
fclose(fid);

% %store transtion records in files
% file1=[sId,'_transToRotation.txt'];
% fid=fopen(file1,'a+');
% for k = 1:size(exactRotationTime,1)
%     for r=1:size(exactRotationTime(k,:),2)
%         fprintf(fid,'%g',exactRotationTime(k,r));
%         fprintf(fid,'\t');
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);
% file2=[sId,'_transToStationary.txt'];
% fid=fopen(file2,'a+');
% for k = 1:size(exactQuitTime,1)
%     for r=1:size(exactQuitTime(k,:),2)
%         fprintf(fid,'%g',exactQuitTime(k,r));
%         fprintf(fid,'\t');
%     end
%     fprintf(fid,'\n');
% end
% fclose(myfile1);
% fclose(myfile2);


%############################
%Clear the screen.
% Eyelink('StopRecording');
% Eyelink('CloseFile');
% try
%     fprintf('Receiving data file ''%s''\n', [sId '.edf'] );
%     status=Eyelink('ReceiveFile', edfFile, pwd, 1);
%
%     if status > 0
%         fprintf('ReceiveFile status %d\n', status);
%     end
%     if 2==exist(edfFile, 'file')
%         fprintf('Data file ''%s'' can be found in ''%s''\n', [sId '.edf'], pwd );
%     end
% catch rdf
%     fprintf('Problem receiving data file ''%s''\n', [sId '.edf'] );
%     rdf;
% end
% Eyelink('Shutdown');
%******************* Eyelink: end
%############################

sca;



