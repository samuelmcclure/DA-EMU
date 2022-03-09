% Three stimulus visual oddball task========================================
% Instruction: Press Space Bar as soon as possible for the target stimulus
% target stimulus: Block 1 = Blue, Block 2 = Red, Block 3 = Green
% Press DownArrow key to start the next block of trials.
% To force quit the experiment: ctrl + c -> command + 0 -> type 'sca'+Enter
% To temporarily stop the experiment, Press and hold Esc key until 'Please
% wait' message on the screen. Press Esc key again to resume the experiment
% 2022-03-04
% GYB
%
%
% Psychtoolbox does not interact with Windows 10 well when drawing fonts.
% This seems to be due to the fact that Windows deals with higher
% resolution displays by scaling fonts to a larger size. Within
% Psychtoolbox this causes fonts to be larger than they ought to be. In
% some text drawing functions, the text is clipped at the top.
% One solution is to turn off font scaling. This can be done by adjusting
% the setting on Windows in the Settings app:
%     Settings > System > Display > Scale and layout --> 100%
% You will probably also want to adjust the screen resolution because
% turning off scaling on a high resolution monitor means that the
% application fonts become miniscule.
% 2022-03-08
% SMM
%==========================================================================

%% THREESTIMVISUALODDBALL function
function [] = ThreeStimVisualOddball_Win

global session;
% InitSession sets a number of variables that define how the task is
% displayed and run. Ideally, you should only need to change the code in
% InitSession to make any necessary tweaks. The code in this file should
% remain fixed across computers.
session = InitSession();
KbName('UnifyKeyNames');

session.subject_id = input('Subject number (please use digits):   ', 's');
session.subject_gender = input('1:Woman, 2:Man, 3:No report:   ', 's');
session.logfile_basename = strcat( ...
        'ThreeStimOddBall_', session.subject_id, '_', ...
        datestr(now, 'mm_dd_HH_MM'));
% make sure that the log file directory exists
mkdir(pwd, 'log');

% Create a connection to the GazePoint eye tracker if the gazepoint should
% be enabled
if session.gazepoint_enable,
    addpath(session.gazepoint_path);

    % this opens a separate Matlab window that maintains a continuous
    % connection to Gazepoint
    session.session1_client = ConnectToGP3;

    % we use the same root filename for task and eyetracking data. the
    % gazepoint data go in at .txt file with the ending 'gazepoint'
    session.gpFileName = strcat(session.logfile_basename, '_gazepoint.txt');
    session.gpFileName = fullfile(pwd, 'log', session.gpFileName);
    ExecuteRecordGP3Data(session.session1_client, session.gpFileName);
end;

% define trials
% NOTE: THIS CODE SHOULD EVENTUALLY BE FIXED. THERE IS NO GUARANTEE THAT
% SESSION.TRIALS_PER_BLOCK EQUAL THE SUM OF NUM_STANDARD + NUM_DEVIANT +
% NUM_TARGET. THERE IS THE POSSIBILITY THAT ROUNDING MAY CAUSE THE NUMBER
% OF TRIALS TO BE OFF BY ONE. 
% I TRIED TO AVOID POSSIBLE ERRORS BY USING NUMEL(SESSION.ALL_TRIALS)
% INSTEAD OF SESSION.TRIALS_PER_BLOCK BUT THIS IS UGLY AND SHOULD BE FIXED.
num_standard = round(session.trials_per_block * session.pct_standard);
num_deviant = round( (session.trials_per_block-num_standard)/2 );
num_target = num_deviant;
standard_trials = ones(1, num_standard);
deviant_trials = 2 * ones(1, num_deviant);
target_trials = 3 * ones(1, num_target);
session.all_trials = [standard_trials, deviant_trials, target_trials];

RunTesting();
% END OF THREESTIMVISUALODDBALL_WIN FUNCTION

%% RUNTESTING function
function [] = RunTesting
global session
    
% initialize data
trial_data = zeros(numel(session.all_trials), 6);  

session.output = ones(session.num_trials, 13);
session.output = -1*session.output;
subject_id = str2double(session.subject_id);
session.output(:,1) = subject_id;

subject_gender = str2double(session.subject_gender);
session.output(:,3) = subject_gender;

%Set one Block
trial = 0;

for tp = 1:length(session.all_trials),
    n_object = 1;
    
    trial_data(trial+1,1) = trial+1; %trial number
    trial_data(trial+1,2) = n_object; % number of objects
    trial_data(trial+1,3) = session.all_trials(tp); % which bin
   
    trial = trial + 1;  
end;

% SHOW INSTRUCTIONS=======================================================
%Screen('Preference', 'SkipSyncTests', 1);
[w, rect] = Screen('OpenWindow', session.monitor, session.background_color);

% Keep track of where the center of the screen is. This is used to place
% the stimuli and fixation point. Having a variable named 'c' is awful, but
% it is used often and is kept short to maintain legible code.
c.x = rect(3)/2;
c.y = rect(4)/2;

% for convenience, remember the rect to draw the stimuli and fixation
% points
r = session.target_size;
stim_rect = [c.x-r, c.y-r, c.x+r, c.y+r];
r = session.fixation_size;
fixation_rect = [c.x-r, c.y-r, c.x+r, c.y+r];    

Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% define how the text is displayed
if ~IsLinux,
    Screen('TextFont', w, session.font_name);
    Screen('TextSize', w, session.font_size);
    Screen('TextStyle', w, session.font_style); % 1+2 is bold and italic
end;
    
%Trial Begins==============================================================
trial_number = 0;
HideCursor();

for block = 1:session.num_blocks,
    if session.gazepoint_enable,
        SendMsgToGP3(session.session1_client,['block_start,' num2str(block)]);
    end;

    %Instruction
    DrawFormattedText(w, ['Target = ' session.target_name{block}], ...
        'center', c.y-session.target_size*1.25, [0 0 0]);
    Screen('FillOval', w, session.target_color{block}, stim_rect); %stimulus
    Screen('FillOval', w, session.fixation_color, fixation_rect);   
    Screen('Flip',w);
    KbWait(-1);

    trial_order = randperm( numel(session.all_trials) );

    for trial = 1 : numel(session.all_trials)
        trial_index = trial_order(trial);
        trial_number = trial_number + 1;
        % standard/deviant/target
        condition = trial_data(trial_index, 3); % trial type: standard, deviant, target
        if condition == 1, % standard
            stim_color = session.standard_color{block};
        elseif condition == 2, % deviant
            stim_color = session.deviant_color{block};
        else, % target
            stim_color = session.target_color{block};
        end;

        if trial == 1,
            % fixation
            Screen('FillOval', w, session.fixation_color, fixation_rect); % fixation dot    
            Screen('Flip',w);
            WaitSecs(session.ITI); % 1s  
        end;

        % Sample Period
        if session.gazepoint_enable,
            SendMsgToGP3(session.session1_client,['trial_start,' num2str(trial) ',' num2str(condition)]);
        end;
        Screen('FillOval', w, stim_color, stim_rect); %stimulus    
        Screen('Flip',w);
        time_stimulus_onset = GetSecs;
        WaitSecs(session.sample);
    
        % Fixation
        Screen('FillOval', w, session.fixation_color, fixation_rect); % fixation dot      
        Screen('Flip',w);
        response_window_onset = GetSecs;
        [secs, keyCode, ~] = KbWait([], 0, response_window_onset+session.blank);
        
        if session.gazepoint_enable,
            keys = '';
            for k = find(keyCode),
                keys = [keys ',' num2str(k)];
            end;
            SendMsgToGP3(session.session1_client,['trial_response' keys]);
        end;
        time_after_response = GetSecs;
    
        % add remaining time of blank sreen
        response_time = time_after_response - response_window_onset;
        if response_time <= session.blank,
            Screen('FillOval', w, session.fixation_color, fixation_rect);       
            Screen('Flip',w);
            % wait out the rest of the response time window
            WaitSecs(session.blank - response_time);
        end;

        if keyCode( KbName('space') ), % space bar
            RT = time_after_response - time_stimulus_onset;
            if condition == 3,
                ACC = 1; % Hit
            else,
                ACC = 0; % False Alarm   
            end;

        elseif keyCode( KbName('ESCAPE') ), % escape key
            resume = false;
            while ~resume,
                DrawFormattedText(w, 'Please wait...', ...
                    'center', c.y-session.target_size*1.25, [0 0 0]);
                Screen('Flip',w);
                [keyIsDown, ~, keyCode] = KbCheck(-1);
                % resume when the escape key is pressed a second time
                resume = keyCode(KbName('ESCAPE'));
                
                WaitSecs(0.1);
            end; % end of trial resume

            % skip saving the trial info
            continue;
        else,
            RT = 0;
            if condition == 3,
                ACC = 0; % Miss
            else,
                ACC = 1; % Correct Rejection   
            end;
        end;
      
        session.output(trial_number, 4) = condition; % 
        session.output(trial_number, 5) = ACC;
        session.output(trial_number, 6) = RT; % 
        session.output(trial_number, 7) = block; % 
        
        SaveMatrix();
    end; %End of one block
    
    %Short delay at the end of each block
    %Fixation
    Screen('FillOval', w, session.fixation_color, fixation_rect);    
    Screen('Flip', w);
    WaitSecs(session.ITI); %1s

    if session.gazepoint_enable,
        SendMsgToGP3(session.session1_client,['block_end,' num2str(block)]);
    end;

    % progress report
    resume_block = false;
    while ~resume_block,
        % Give progress report. 
        progress_string = num2str( round((block/session.num_blocks)*100) );
        progress_string = strcat(progress_string, '% done.');
        
        if block == session.num_blocks,
            progress_string = strcat(progress_string, '\n\nThank you!')
        else,
            progress_string = strcat(progress_string, '\n\nTake a short break.');       
        end;
        DrawFormattedText(w, progress_string, 'center', 'center', [0 0 0]);
        Screen('Flip',w);

        [keyIsDown, ~, keyCode] = KbCheck(-1);
        resume_block = keyCode( KbName('DownArrow') );
    end;
    WaitSecs(session.ITI);
end; %End of all trials

if session.gazepoint_enable,
    %% Stop collecting data in client2
    fprintf('Stop recording\n')
    SendMsgToGP3(session.session1_client,'STOP_EYETRACKER');

    %% Clean-up socket
    CleanUpSocket(session.session1_client);
    % This closes all open files. Note that this could cause problems, but
    % it is unclear how to close the output file for Gazepoint recording 
    % selectively.
    fclose all; 
end;

WaitSecs(0.5);
Screen('Close', w); 
ShowCursor;

% END RUNTESTING FUNCTION


%% SAVEMATRIX FUNCTION
function [] = SaveMatrix()

global session;
subject_id = num2str(session.subject_id);
fn = strcat(session.logfile_basename, '.xls');
fn = fullfile(pwd, 'log', fn);

for p=1:session.num_trials,
   session.output(p,2) = p; 
end;
       
fid = fopen(fn, 'w');
fprintf(fid, '%s\t %s\t %s\t %s\t %s\t %s\t  %s\t \n', ...
    'Subject_Num',...
    'Trial_Num',...
    'Gender',...
    'Condition',...
    'ACC',...
    'RT',...
    'Block');
     
for p=1:session.num_trials
    fprintf(fid, '%4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t \n', ...
        session.output(p, 1),...
        session.output(p, 2),...
        session.output(p, 3),...
        session.output(p, 4),...
        session.output(p, 5),...
        session.output(p, 6),...
        session.output(p, 7));
end;
fclose(fid);

% END SAVEMATRIX FUNCTION
