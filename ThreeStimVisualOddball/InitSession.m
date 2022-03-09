function session = InitSession()

%% Number of trials, blocks, etc.
session.trials_per_block = 60;
session.pct_standard = 0.7;

session.standard_color = { [255,0,0], [0,200,0], [0,0,255] };
session.deviant_color = { [0,200,0], [0,0,255], [255,0,0] };
session.target_color = { [0,0,255], [255,0,0], [0,200,0] };
session.target_name = { 'Blue', 'Red', 'Green' };

session.num_blocks = numel(session.standard_color);
session.num_trials = session.num_blocks * session.trials_per_block;

%% Display preferences (colors, sizes, etc.)
% monitor number
session.monitor = 0;

% target size in pixels (radius)
session.target_size = 100;

% size of fixation circle in pixels (radius)
session.fixation_size = 2;

% font
session.font_name = 'Arial';
session.font_size = 50;
session.font_style = 1+2; % bold and italic: see http://psychtoolbox.org/docs/Screen-TextStyle

% colors
session.background_color = [125 125 125];
session.fixation_color = [0 0 0];

%durations
session.sample = 0.3; %1
session.blank = 0.7; %0.5
session.ITI = 1;

%% Gazepoint setup
% You must have the gazepoint-matlab-toolbox stored on your computer to use
% Gazepoint to track eye movements within this task. You can clone the 
% reposity from this address:
% https://github.com/RingoHHuang/gazepoint-matlab-toolbox
% The variable session.gazepoint_path should match the location of your
% cloned repository.
session.gazepoint_enable = false;
session.gazepoint_path = 'C:\Users\McClureLab\Desktop\ASU-Banner\gazepoint-matlab-toolbox\GP3_Functions\';