%Three stimulus visual oddball task========================================
%Instruction: Press Space Bar as soon as possible for the target stimulus
% target stimulus: Block 1 = Blue, Block 2 = Red, Block 3 = Green
% Press DownArrow key to start the next block of trials.
%To force quit the experiment: ctrl + c -> command + 0 -> type 'sca'+Enter
%To temporarily stop the experiment, Press and hold Esc key until 'Please
%wait' message on the screen. Press Esc key again to resume the experiment
%2022-03-04
%GYB
%==========================================================================

function [] = ThreeStimVisualOddball

global subject_id
global sample
global blank
global num_trials
global gray
global white
global black
global ITI
global nBlock
global trialPerBlock
global trialType
global tarSz
global subject_gender

subject_id = input('Subject number (please use digits):   ', 's');
subject_gender = input('1:Woman, 2:Man, 3:No report:   ', 's');

Ratios = input('Standard (e.g., 0.7):   ', 's');
Ratios = str2double(Ratios);

% target size in pixel (radius)
tarSz = 100;

% basic color
gray = [125 125 125];
black = [0 0 0];
white = [225,225,225];

%durations
sample = 0.3; %1
blank = 0.7; %0.5
ITI = 1;

%number of trials
nBlock = 3;
trialPerBlock = 60;
num_trials  = trialPerBlock * nBlock;
standard = ones(1,trialPerBlock * Ratios);
dev_1 = 2*ones(1,round(trialPerBlock*(1-Ratios)/2));
dev_2 = 3*ones(1,round(trialPerBlock*(1-Ratios)/2));

trialType = [standard,dev_1,dev_2];

runTesting()
end

function [] = runTesting
global subject_id subject_gender num_trials rect  gray black  sample trialType   nBlock ...
    trialPerBlock blank   tarSz  ITI
    
%create Trial data
% usedTrials = zeros(num_trials);
trialData = zeros(trialPerBlock, 6);  

output = ones(num_trials, 13);
output = -1*output;
subject_id = str2double(subject_id);
output(:,1) = subject_id;

subject_gender = str2double(subject_gender);
output(:,3) = subject_gender;


%Set one Block
trial = 0;

  for tp = 1:length(trialType) %4
   n_object = 1;
      trialData(trial+1,1) = trial+1; %trial number
      trialData(trial+1,2) = n_object; % number of objects
      trialData(trial+1,3) = trialType(tp); % which bin
   trial = trial + 1;
  end

% SHOW INSTRUCTIONS=======================================================
curMon= 0;
 %Screen('Preference', 'SkipSyncTests', 1);
%[w,rect]=Screen('OpenWindow',curMon,gray, [0,0,1600,670]);
[w,rect]=Screen('OpenWindow',curMon,gray);

    %target location
    targx = rect(3)/2 ;
    targy = rect(4)/2 ;
    
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
%Text assignment
    if IsLinux==0
        Screen('TextFont',w, 'Arial');
        Screen('TextSize',w, 50);
        Screen('TextStyle', w, 1+2);
    end
    

%Trial Begins==============================================================
trialNo = 0;
%firstTrial = 0;
HideCursor()

for block = 1:nBlock %repeat multiple times
    
    %Instruction
    if block == 1
    Screen('DrawText', w, 'Target = Blue ', rect(3)/2-150, rect(4)/2-200,[0 0 0]);    
    Screen('FillOval', w, [0,0,255],[targx-tarSz,targy-tarSz,targx+tarSz,targy+tarSz]);%stimulus    
    elseif block == 2
    Screen('DrawText', w, 'Target = Red ', rect(3)/2-150, rect(4)/2-200,[0 0 0]);    
    Screen('FillOval', w, [255,0,0],[targx-tarSz,targy-tarSz,targx+tarSz,targy+tarSz]);%stimulus    
    elseif block == 3
    Screen('DrawText', w, 'Target = Green ', rect(3)/2-150, rect(4)/2-200,[0 0 0]);    
    Screen('FillOval', w, [0,255,0],[targx-tarSz,targy-tarSz,targx+tarSz,targy+tarSz]);%stimulus    
    end
    Screen('FillOval', w, black,[rect(3)/2-2,rect(4)/2-2,rect(3)/2+2,rect(4)/2+2]);%SMallDot    
    Screen('Flip',w);
    %KbWait(-1)
    KbWait;

trialorder = randperm(trialPerBlock);

for trial = 1 : trialPerBlock
 
    triIndex = trialorder(trial);
    trialNo = trialNo + 1;
    
    
 % standard/deviant/target
 condition = trialData(triIndex,3); % trial type: standard, deviant, target
 
     if condition == 1 % standard 
        if block == 1
         stimCol = [255,0,0]; % red
        elseif block == 2
         stimCol = [0,200,0]; % green  
        elseif block == 3
         stimCol = [0,0,255]; % blue 
        end
     elseif condition == 2 % deviant
        if block == 1
         stimCol = [0,200,0]; %green
        elseif block == 2
         stimCol = [0,0,255]; %blue 
        elseif block == 3
         stimCol = [255,0,0]; %red 
        end
     else % target
        if block == 1
         stimCol = [0,0,255]; % blue
        elseif block == 2
         stimCol = [255,0,0]; % red  
        elseif block == 3
         stimCol = [0,200,0]; % green 
        end
     end

    if trial == 1
    %Fixation
    Screen('FillOval', w, black,[rect(3)/2-2,rect(4)/2-2,rect(3)/2+2,rect(4)/2+2]);%SMallDot    
    Screen('Flip',w);
    WaitSecs(ITI); %1s
    end

 %Sample Period
    Screen('FillOval', w, stimCol,[targx-tarSz,targy-tarSz,targx+tarSz,targy+tarSz]);%stimulus    
    Screen('Flip',w);
    Onset = GetSecs;
    WaitSecs(sample)
    
    Screen('FillOval', w, black,[rect(3)/2-2,rect(4)/2-2,rect(3)/2+2,rect(4)/2+2]);%SMallDot       
    Screen('Flip',w);
    curTime = GetSecs;
    %[secs, keyCode,~ ] = KbWait(-1,0,curTime+blank);
    [secs, keyCode,~ ] = KbWait([],0,curTime+blank);
    PostResp = GetSecs;
    
    % add remaining time of blank sreen
    TimeDiff = PostResp - curTime;
    if TimeDiff <=blank % 0.7 blank
    Screen('FillOval', w, black,[rect(3)/2-2,rect(4)/2-2,rect(3)/2+2,rect(4)/2+2]);%SMallDot       
    Screen('Flip',w);
    WaitSecs(blank-TimeDiff);
    end

            if keyCode(44) % space bar
            RT = PostResp - Onset;
                if condition == 3
                ACC = 1; % Hit
                else
                ACC = 0; % False Alarm   
                end

            elseif keyCode(41) % escape key
                resume=0;
                while resume==0
                    Screen('DrawText', w, 'Please wait....', rect(3)/2-150, rect(4)/2-200,[0 0 0]);    
                    Screen('Flip',w);
                    [ keyIsDown,~, keyCode ] = KbCheck(-1);
                    if keyCode(41) % escape key
                    resume=1;
                    end
                    WaitSecs(0.1);
                end % end of trial resume
            else
                RT = 0;
                if condition == 3
                ACC = 0; % Miss
                else
                ACC = 1; % Correct Rejection   
                end

            end
      
     output(trialNo, 4) = condition; % 
     output(trialNo, 5) = ACC;
     output(trialNo, 6) = RT; % 
     output(trialNo, 7) = block; % 
 
     saveMatrix();

end %End of one block
    %Short delay at the end of each block
    %Fixation
    Screen('FillOval', w, black,[rect(3)/2-2,rect(4)/2-2,rect(3)/2+2,rect(4)/2+2]);%SMallDot    
    Screen('Flip',w);
    WaitSecs(ITI); %1s
    
    % progress report
    Blockresume =0;
    while Blockresume ==0
    % Give progress report. 
    Progress = num2str(round((block/nBlock)*100));
    DispProgress = strcat(Progress,' % is done.');
    Screen('DrawText', w, DispProgress, rect(3)/2-300, rect(4)/2-100,[0 0 0]);    
    if block == nBlock
    Screen('DrawText', w, 'Thank you.', rect(3)/2-300, rect(4)/2,[0 0 0]);    
    else
    Screen('DrawText', w, 'Take a short break.', rect(3)/2-300, rect(4)/2,[0 0 0]);        
    end
    Screen('Flip',w);
        [ keyIsDown,~, keyCode ] = KbCheck(-1);
        if keyCode(81) % down arrow key
        Blockresume=1;
        end
    end
end %End of whole trials
WaitSecs(0.5); 

function [] = saveMatrix()
        subject_id = num2str(subject_id);
        fn = strcat('ThreeStimOddBall_', subject_id, '_',datestr(now, 'mm_dd'),'_', '.xls');
        
        for p=1:num_trials
           output(p,2) = p; 
        end
       
        fid = fopen(fn, 'w');
        fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t  %s\t \n', ...
        'Subject_Num',...
        'Trial_Num',...
        'Gender',...
        'Condiiton',...
        'ACC',...
        'RT',...
        'Block');
     
        for p=1:num_trials
            fprintf(fid, '%4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t \n', ...
            output(p, 1),...
            output(p, 2),...
            output(p, 3),...
            output(p, 4),...
            output(p, 5),...
            output(p, 6),...
            output(p, 7));
        end
        fclose(fid);
end

    Screen('Close',w) 
    ShowCursor
end