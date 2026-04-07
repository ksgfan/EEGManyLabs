%%

p.backColor = [128, 128, 128];
p.white = [200 200 200];
p.gray = [120 120 120];

prefs = p;

colorBackground = prefs.backColor; %[128 128 128];

minStimDuration = 0.095; %Minimum time to display stim before response is allowed (saccades take ~200ms)
minFixDuration = 0.195; % Time after spacebar indicating fixation before target disappears
isiDuration = 1.2;
jitter = 0.4;

colorProbe = [255 0 0]; % color of the saccade target

fixationSize = 30; %angle2pix(prefs, 0.6); %Size of the outermost ring of the fixation bullseye in angles
fixationStep = (fixationSize - 0.5) / 3; % Step size of the bullseye (fixationSize - 3*fixationStep = fixationStep)

%
screenX = 800; % Screen width
screenY = 600; % Screen height
cx = screenX / 2; % Center X
cy = screenY / 2; % Center Y

% Set up the figure
figure;
hold on;
axis equal;
set(gca, 'Color', [colorBackground/255]); % Set the background color

% Plot the fixation dots
colors = {[1 1 1], [0 0 0], [1 1 1], [0 0 0]}; % Alternating white and black
sizes = [fixationSize, fixationSize-fixationStep, fixationSize-(2*fixationStep), fixationSize-(3*fixationStep)];

for i = 1:4
    scatter(cx, cy, sizes(i)^2, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none');
end

% Plot the saccade target (example position)
currentLoc = [cx + 100, cy]; % Example target location

% Plot the saccade target dots
for i = 1:4
    scatter(currentLoc(1), currentLoc(2), sizes(i)^2, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none');
end

% Set the figure limits
xlim([0 screenX]);
ylim([0 screenY]);
set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
hold off;


%% video

% Video parameters
videoFileName = 'saccade_stimulus.mp4';
videoFrameRate = 1; % Frames per second
nFrames = 3; % Total number of frames (adjust as needed)
pauseDurationInFrames = 1
% Create a video writer object with 'MPEG-4' format
v = VideoWriter(videoFileName, 'MPEG-4');
v.FrameRate = videoFrameRate;
open(v);

% Set up figure for recording
figure;
hold on;
axis equal;
set(gca, 'Color', [colorBackground/255]); % Set the background color
xlim([0 screenX]);
ylim([0 screenY]);
set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes

% Alternating saccade target locations
targetOffsets = [-150, 150]; % Left and right offsets from the center
targetOffsets2 = [-300, 300]; % Left and right offsets from the center
frameStep = 2; % Frame interval for switching target position

for frame = 1:nFrames
    % Clear the figure
    clf;
    hold on;
    set(gca, 'Color', [colorBackground/255]); % Set the background color
    axis equal;
    
        axis equal;
    set(gca, 'Color', [colorBackground/255]); % Set the background color
    xlim([0 screenX]);
    ylim([0 screenY]);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    
    % Plot the fixation dots (center)
    for i = 1:4
        scatter(cx, cy, sizes(i)^2, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none');
    end
    
    axis equal;
    set(gca, 'Color', [colorBackground/255]); % Set the background color
    xlim([0 screenX]);
    ylim([0 screenY]);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    
    pause(1)
    frameCapture = getframe(gcf);
    for repeat = 1:1% pauseDurationInFrames
        writeVideo(v, frameCapture);
    end
    
    randomNumber = randi([0, 1]);
    
    if randomNumber == 0
        % Determine the target location based on the frame number
        currentOffset = targetOffsets(mod(floor(frame / frameStep), 2) + 1);
    else
        currentOffset = targetOffsets2(mod(floor(frame / frameStep), 2) + 1);
    end
        
    currentLoc = [cx + currentOffset, cy];
    
    % Plot the saccade target dots
    for i = 1:4
        scatter(currentLoc(1), currentLoc(2), sizes(i)^2, 'MarkerFaceColor', colors{i}, 'MarkerEdgeColor', 'none');
    end
    
    axis equal;
    set(gca, 'Color', [colorBackground/255]); % Set the background color
    xlim([0 screenX]);
    ylim([0 screenY]);
    set(gca, 'XColor', 'none', 'YColor', 'none'); % Hide axes
    
    pause(2)
    for repeat = 1:pauseDurationInFrames
        writeVideo(v, frameCapture);
    end
    
    % Capture the frame
    frameCapture = getframe(gcf);
    writeVideo(v, frameCapture);
    
    
    hold off;
end

close(v);




%%
