%Color wheel memory task with set size 2===================================
%Instruction: Remember two colors and click on the colorwheel using a mouse
%to report the color of the probed object.
%Participants will receive +5 points if response error for a given trial is
%less than 15 degrees.
%This code requires ASUcolor.csv file.
%To force quite the experiment: ctrl + c -> command + 0 -> type 'sca'+Enter
%2022-03-04
%GYB
%==========================================================================
function [] = ColorWheelMemory_Set2_V1_Mac %main function

global subject_id
global rectx
global recty
global sample
global blank
global num_trials
global RadiusPosition
global position
global ColorRGBs
global OuterRadius
global InnerRadius
global width
global NColors
global gray
global SetSize
global nBlock
global ColorSubsetIDs
global NColorSamples
global TrialPerBlock
global gender
global ColDiff
global CWCCW

subject_id = input('Subject number (numeric):   ', 's');
gender = input('1: Woman, 2: Man, 3: No response: ','s');
ColorSet  = csvread('ASUcolor.csv',1,0); % 180 RGBs from CIELAB space

NColors=180; % total number of RGBs in ColorSet
MinColDiff = 12; % 1 color unit = 2 degree on the color wheel
ColorSubsetIDs  = 1:MinColDiff:NColors;  % colors to be used as targets
NColorSamples = length(ColorSubsetIDs);
ColorRGBs = round(ColorSet(1:NColors,:));

%DrawColorWheel
OuterRadius = 275;
InnerRadius = 220;
width = 5;

%square half size
rectx = 32; %size of the object in pixels
recty = 32; %size of the object in pixels
RadiusPosition = 148; % radius of an invisible circle on which the objects willl be presented.
position = 8; % 8 equally spaced locations for the objects

gray = [50 50 50];

%SetSize
SetSize = 2;

% color similarity: 1 = 2 degree on the wheel

ColDiff = [15, 45]; % small differene = 30 degree. large difference = 45 degree
CWCCW = [-1,1];
%durations
sample = 0.2;
blank = 0.9;

%number of trials
HowManySet = length(SetSize);
nBlock = 3;
Nsim = 2; % Non-target color similarity: similar vs dissimilar
NCW = 2; % Non-target color relative to the target: clockwise or counter-clockwsie from the target

TrialPerBlock =NColorSamples *NCW* Nsim* HowManySet; 
num_trials =  TrialPerBlock * nBlock;


runTesting()
end

function [] = runTesting
global subject_id rectx recty sample blank num_trials rect RadiusPosition position ColorRGBs NColorSamples ...
       OuterRadius InnerRadius width NColors gray SetSize nBlock ColorSubsetIDs TrialPerBlock gender ColDiff CWCCW


%create Trial data
trialData = zeros(TrialPerBlock, 5);  

output = ones(num_trials, 18);
output = -1*output;
subject_id = str2double(subject_id);
output(:,1) = subject_id;

gender = str2double(gender);
output(:,17) = gender;

n_object = SetSize(1);

ClickPosition = zeros(num_trials,2); % x-y coordiates of mouse click
UsedColorIndex = zeros(num_trials,n_object);
UsedColorRGB = -1*ones(num_trials,n_object*3);
    
trial = 0;
    for q = 1 : NColorSamples
        for  sim = 1:2 % similar vs dissimilar
            for relPos = 1:2 % clockwise vs counter cw
                
                  trialData(trial+1,1) = trial+1; %trial number
                  trialData(trial+1,2) = n_object; % number of objects
                  trialData(trial+1,3) = ColDiff(sim); % color value of the first targt
                  trialData(trial+1,4) = ColorSubsetIDs(q); % color value of the first targt
                  trialData(trial+1,5) = CWCCW(relPos); % color value of the first targt
                  
               trial = trial + 1;
            end
        end
    end
  
HideCursor;
% SHOW INSTRUCTIONS=======================================================
whichScreen = 0;
[w,rect]=Screen('OpenWindow',whichScreen,gray);
%Text assignment
    if IsLinux==0
        Screen('TextFont',w, 'Courier');
        Screen('TextSize',w, 50);
        Screen('TextStyle', w, 1+2);
    end
    
DrawFormattedText(w,'Main session - color memory set 2','center', 'center',[255 255 255]);
DrawFormattedText( w,'Press SpaceBar to begin','center', rect(4)/2+100,[255 255 255]);


Screen('Flip',w); 
KbWait;
WaitSecs(0.3);


%Trial Begins==============================================================
trialNo = 0;

ColorPalette = -1*ones(num_trials,n_object); %This is color used in a trial from the full color set

TotalScore = 0; % accumulate score

for block = 1:nBlock
    
    trialorder = randperm(TrialPerBlock);
    
for trial = 1 : TrialPerBlock

    triIndex = trialorder(trial);
    trialNo = trialNo + 1;

    % stimulus color(RGBs) used in a given trial
    color = zeros(n_object,3);
    
    %Pick random color value from new color index matrix.
    TargetColor = trialData(triIndex,4); % target color index
    
    % similar or dissimilar
    colDiffcond = trialData(triIndex,3);
    whichSide = trialData(triIndex,5);
        
    NTColor = TargetColor + whichSide*colDiffcond;

    NTColor = wrap(NTColor,1,181); % wrap color ID: colors should be in [1,180]
    
    colSet = [TargetColor,NTColor]; % stimulus color IDs
        
    for i = 1 : n_object            
       color(i,:) = ColorRGBs(colSet(i),:); % This is the color used in a trial, This ColorWheel is newly ordered colorwheel
       ColorPalette(trialNo,i) = colSet(i); %Absoluet index of color
    end
        
    % stimulus location    
    x = zeros(n_object, 1);
    y = zeros(n_object, 1);
    StimPosition = NRandPerm(position, n_object);
    RecordPosition =zeros(6,1);
     
    for i = 1 : n_object %1 is taret position
    x(i) = rect(3)/2+RadiusPosition*cos(StimPosition(i)*2*pi/position);
    y(i) = rect(4)/2+RadiusPosition*sin(StimPosition(i)*2*pi/position);
    RecordPosition(i,1) = StimPosition(i);
    end

 %Fixation Period       
    %fixation
    Screen('DrawLine', w, [255 255 255], (rect(3)/2)-10, rect(4)/2, rect(3)/2+10, rect(4)/2,5);
    Screen('DrawLine', w, [255 255 255], rect(3)/2, (rect(4)/2)-10, rect(3)/2, (rect(4)/2)+10,5);    
    Screen('Flip',w);
    WaitSecs(0.5);
    
 %Sample Period
    %fixation
    Screen('DrawLine', w, [255 255 255], (rect(3)/2)-10, rect(4)/2, rect(3)/2+10, rect(4)/2,5);
    Screen('DrawLine', w, [255 255 255], rect(3)/2, (rect(4)/2)-10, rect(3)/2, (rect(4)/2)+10,5);    
    
    add = 0;
    for q = 1:n_object
        Screen('FillRect',w,color(q,:),[x(q)-rectx,y(q)-recty,x(q)+rectx,y(q)+recty]); 
        UsedColorIndex(trialNo,q) = ColorPalette(trialNo,q);
        UsedColorRGB(trialNo,q+add) = color(q,1); % R
        UsedColorRGB(trialNo,q+add+1) = color(q,2); % G
        UsedColorRGB(trialNo,q+add+2) = color(q,3); %B
        add = add+2;
    end
    Screen('Flip',w);
    WaitSecs(sample);
  
 %Delay period
   %fixation
   Screen('DrawLine', w, [255 255 255], (rect(3)/2)-10, rect(4)/2, rect(3)/2+10, rect(4)/2,5);
   Screen('DrawLine', w, [255 255 255], rect(3)/2, (rect(4)/2)-10, rect(3)/2, (rect(4)/2)+10,5);    
   Screen('Flip',w);
   WaitSecs(blank);
    
   
 %Test Period
ColorAngles = -1*ones(1,n_object);
SetMouse(rect(3)/2,rect(4)/2,whichScreen)    %initial location of the cursor    
clicks = 0;
    % Initial position of cursor
    ShowCursor(whichScreen)
    % set initial value 
    st =NRandPerm(NColors,1);

    while clicks ==0

    neworder = wrap(st:(st+179),1,181);
    ColorWheel = ColorRGBs(neworder,:);
    
    %wheel
    for i = 1: NColors
       Screen('DrawLine', w, ColorWheel(i,:),rect(3)/2, rect(4)/2,...
       rect(3)/2+OuterRadius*cos(i*2*pi/NColors),rect(4)/2+OuterRadius*sin(i*2*pi/NColors),width);
       Screen('DrawLine', w, ColorWheel(i,:),rect(3)/2, rect(4)/2,...
       rect(3)/2+OuterRadius*cos((i+0.5)*2*pi/NColors),rect(4)/2+OuterRadius*sin((i+0.3)*2*pi/NColors),width);
        Screen('DrawLine', w, ColorWheel(i,:),rect(3)/2, rect(4)/2,...
       rect(3)/2+OuterRadius*cos((i+0.6)*2*pi/NColors),rect(4)/2+OuterRadius*sin((i+0.6)*2*pi/NColors),width);
    end
    Screen('FillOval', w, gray, [rect(3)/2-InnerRadius rect(4)/2-InnerRadius rect(3)/2+InnerRadius rect(4)/2+InnerRadius]);
    Screen('FrameOval', w, gray, [rect(3)/2-(OuterRadius+3) rect(4)/2-(OuterRadius+3) rect(3)/2+(OuterRadius+3) rect(4)/2+(OuterRadius+3)],10);
    
    %objects
    for q = 1:n_object
        Screen('FrameRect',w,[255 255 255],[x(q)-rectx,y(q)-recty,x(q)+rectx,y(q)+recty], 1); 
    end
    %target object
    Screen('FrameRect',w,[255 255 255],[x(1)-rectx,y(1)-recty,x(1)+rectx,y(1)+recty],8); 
    
    %Fixation
    Screen('DrawLine', w, [255 255 255], (rect(3)/2)-10, rect(4)/2, rect(3)/2+10, rect(4)/2,5);
    Screen('DrawLine', w, [255 255 255], rect(3)/2, (rect(4)/2)-10, rect(3)/2, (rect(4)/2)+10,5);    
    Screen('Flip',w);
    tic
    
    % Get Clicked location
    [cursorx,cursory,button] = GetMouse(whichScreen,0);
    
        if button(1) ==1
        ElapsedT=toc; % response time
        ClickPosition(trialNo,1) = cursorx;
        ClickPosition(trialNo,2) = cursory;
        response(1) = cursorx;
        response(2) = cursory;
        clicks=1;
        end
    
    end % End of while loop for each click
%%  

    % response angle
    AngleFromCenterToClick = atand((response(2) - rect(4)/2)/(response(1) - rect(3)/2));
    
    if AngleFromCenterToClick >=0
       if cursory - rect(4)/2 >=0
       else
          AngleFromCenterToClick = AngleFromCenterToClick + 180;
       end
    else
         if cursory - rect(4)/2 >=0
          AngleFromCenterToClick = 180 + AngleFromCenterToClick;
         else
          AngleFromCenterToClick = 360 + AngleFromCenterToClick;
         end
    end
    
    %record response 
    output(trialNo, 7) = AngleFromCenterToClick; %
       cid = wrap(floor(AngleFromCenterToClick/2),1,181);
    if cid ==0
       cid =180;
    else
    end
    RColorID = neworder(cid);
    
    
    % mark response on the wheel
    Screen('drawline', w, [255 255 255],rect(3)/2, rect(4)/2,...
                          rect(3)/2+(OuterRadius+30)*cosd(AngleFromCenterToClick),rect(4)/2+(OuterRadius+30)*sind(AngleFromCenterToClick),width);

    %DrawColorWheel
    for i = 1: NColors
       Screen('DrawLine', w, ColorWheel(i,:),rect(3)/2, rect(4)/2,...
       rect(3)/2+OuterRadius*cos(i*2*pi/NColors),rect(4)/2+OuterRadius*sin(i*2*pi/NColors),width);
       Screen('DrawLine', w, ColorWheel(i,:),rect(3)/2, rect(4)/2,...
       rect(3)/2+OuterRadius*cos((i+0.3)*2*pi/NColors),rect(4)/2+OuterRadius*sin((i+0.3)*2*pi/NColors),width);
        Screen('DrawLine', w, ColorWheel(i,:),rect(3)/2, rect(4)/2,...
       rect(3)/2+OuterRadius*cos((i+0.6)*2*pi/NColors),rect(4)/2+OuterRadius*sin((i+0.6)*2*pi/NColors),width);
    end
    Screen('FillOval', w, gray, [rect(3)/2-InnerRadius rect(4)/2-InnerRadius rect(3)/2+InnerRadius rect(4)/2+InnerRadius]);
    Screen('FrameOval', w, gray, [rect(3)/2-(OuterRadius+3) rect(4)/2-(OuterRadius+3) rect(3)/2+(OuterRadius+3) rect(4)/2+(OuterRadius+3)],10);
    
    %Link between color index(based on non-rotated colorwheel) and rotated color wheel
    ObjectColorNumber = zeros(n_object,1);
    for i = 1 : n_object
    ObjectColorNumber(i,1) = UsedColorIndex(trialNo,i);
    end
    
    % response color
    for obj = 1:n_object
    for j = 1:NColors        
        if neworder(1,j) ==  UsedColorIndex(trialNo,obj)
        ColorAngles(1,obj) = j;
        end
    end
    end

     % Show accuracy feedback    
     currentE =  wrap(ColorAngles(1,1)*2 - AngleFromCenterToClick,-180,180);
       
       if abs(currentE) <= 15
       show = '+5';
       currentScore = 5;
       TotalScore = TotalScore+currentScore;
       else
       show = '+0';    
       currentScore = 0;
       TotalScore = TotalScore+currentScore;
       end
           
      Screen('DrawText', w, show, rect(3)/2-50, rect(4)/2-10,[255 255 255]);      
    
%     % target color position   
%     currentTarget = ColorRGBs(UsedColorIndex(trialNo,1),:); % currently reported color    
%     Screen('DrawLine', w, currentTarget,rect(3)/2, rect(4)/2,...
%        rect(3)/2+(OuterRadius+30)*cos(ColorAngles(1,1)*2*pi/NColors),rect(4)/2+(OuterRadius+30)*sin(ColorAngles(1,1)*2*pi/NColors),width);
%     % non-target color position   
%    currentNT = ColorRGBs(UsedColorIndex(trialNo,2),:); % currently reported color    
%     Screen('DrawLine', w, currentNT,rect(3)/2, rect(4)/2,...
%        rect(3)/2+(OuterRadius+30)*cos(ColorAngles(1,2)*2*pi/NColors),rect(4)/2+(OuterRadius+30)*sin(ColorAngles(1,2)*2*pi/NColors),width);

   
    % compute angular position of target color
    ObjectAngleDegree = zeros(n_object,1);
    for i = 1:n_object
    ObjectAngleDegree(i,1) = convert2deg(ColorAngles(1,i)*2*pi/NColors);
    end
    Screen('Flip',w);
    HideCursor;
    WaitSecs(0.5);
   
    %% Angle of each objects
    output(trialNo, 3) = n_object;
    output(trialNo, 4) = ObjectColorNumber(1,1); %TargetColorPosition
    output(trialNo, 5) = ObjectAngleDegree(1,1); % 
    output(trialNo, 6) = TotalScore; % score
    % output(trialNo, 7) response angle is recorded at line 311
    output(trialNo,8) = RColorID; % report color
    
    for i = 1: n_object
    output(trialNo, i+8) = ObjectAngleDegree(i,1); % D btw center of the object and click                         
    end
    
    output(trialNo, 11) = st; % RotationUnit
    output(trialNo, 12) = ElapsedT; % Responsetime
    output(trialNo, 13) = RecordPosition(1,1);
    output(trialNo, 14) = RecordPosition(2,1);
    output(trialNo, 15) = colDiffcond;
    output(trialNo, 16) = whichSide;
    
    saveMatrix();

end
% ITI
    if block < nBlock
    MaxScore = TrialPerBlock*block*5;
    PercentScore = round(100*TotalScore/MaxScore);
    
    feed = strcat('You earned'," ",num2str(PercentScore),'% of the total score.');
    DrawFormattedText(w,feed,'center', rect(4)/2-100,[255 255 255]);
    DrawFormattedText(w,'Take a short break','center.', rect(4)/2,[255 255 255]);
    DrawFormattedText(w,'Press space bar to begin','center.', (rect(4)/2)+100,[255 255 255]);
    Screen('Flip',w);
    KbWait;
    
    elseif block == nBlock    
    DrawFormattedText(w,'Thank you!','center', 'center',[255 255 255]);
    Screen('Flip',w);
    WaitSecs(0.5)
    
    else % intertrial interval
    Screen('DrawLine', w, [255 255 255], (rect(3)/2)-10, rect(4)/2, rect(3)/2+10, rect(4)/2,5);
    Screen('DrawLine', w, [255 255 255], rect(3)/2, (rect(4)/2)-10, rect(3)/2, (rect(4)/2)+10,5);    
    Screen('Flip',w);
    WaitSecs(1);
    end

end % end of block


function [] = saveMatrix()
        subject_id = num2str(subject_id);
        fn = strcat('ColorWheelMemory_Set2_', subject_id, '_',datestr(now, 'mm_dd'),'_', '.xls');
        
        
        for p=1:num_trials
           output(p,2) = p; 
        end
        

        fid = fopen(fn, 'w');


        fprintf(fid, '%s\t %s\t %s\t %s\t  %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t  \n', ...
        'subject_Num',...
        'trial_Num',...
        'num_objects',...
        'TargetColor',...
        'TargetAngle',...
        'Score',...
        'resp1',...
        'RespCol1',...
        'AngleO1',...
        'AngleO2',...
        'Rotation',...
        'ElapsedT',...
        'O1Location',...
        'O2Location',...
        'Similarity',...
        'NTside',...
        'Gender',...
        'xposition',...
        'yposition',...
        'O1Color',...
        'O2Color');
            
        for p=1:num_trials
    
            fprintf(fid, '%4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t %4d\t \n', ...
            output(p, 1),...
            output(p, 2),...
            output(p, 3),...
            output(p, 4),...
            output(p, 5),...
            output(p, 6),...
            output(p, 7),...
            output(p, 8),...
            output(p, 9),...
            output(p, 10),...
            output(p, 11),...
            output(p, 12),...
            output(p, 13),...
            output(p, 14),...
            output(p, 15),...
            output(p, 16),...
            output(p, 17),...
            ClickPosition(p,1),...
            ClickPosition(p,2),...
            UsedColorIndex(p,1),...
            UsedColorIndex(p,2));
        end

        fclose(fid);

end
    
 
    Screen('Close',w) 
    ShowCursor
end
function wraping = wrap(x, xmin, xmax)
wraping = x - floor((x-xmin)/(xmax-xmin))*(xmax-xmin);
end