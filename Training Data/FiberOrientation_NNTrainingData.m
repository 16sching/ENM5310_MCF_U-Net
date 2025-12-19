%% Final Adapted MATLAB Script: 2D Fibril Angle Analysis with Visualization

clc; clear; close all;

% Get all TIFF files in the current directory
imageFiles = dir('*.tiff');
if isempty(imageFiles)
    error('No TIFF files found in the current directory.');
end

% Image resolution
pixelpernm = 1/21; % 1 pixel = 21 nm

% Initialize a cell array to store results for CSV output
resultsData = {};
resultRowIndex = 1;

%% Define Visualization Colormap
startColor = [1, 1, 1];          % White
endColor = [64, 1, 84] / 255;    % Dark purple
nColors = 256;
customMap = [linspace(startColor(1), endColor(1), nColors)', ...
             linspace(startColor(2), endColor(2), nColors)', ...
             linspace(startColor(3), endColor(3), nColors)'];

%% Step 1: Define Sub-Window Parameters
subwindowSize = [100, 100]; % Size of each sub-window in x,y directions
overlap = 0.5;              % Overlap by 50%
stepFactor = (1 - overlap); 

% Line rotation sweep angles
angles = -90:0.5:90; % Rotation angles of profile line

% Pre-calculate rotation line points (based on sub-window size)
half_size = subwindowSize(1)/2;
x1_rot = half_size - (half_size) * cosd(angles'-90);
y1_rot = half_size - (half_size) * sind(angles'-90);
x2_rot = half_size + (half_size) * cosd(angles'-90);
y2_rot = half_size + (half_size) * sind(angles'-90);
Rot_points = [x1_rot, y1_rot, x2_rot, y2_rot];


%% Step 2, 3, 4, 5: Process Each 2D Image (First Pass, Second Pass, Store Data)
disp('Starting 2D image processing and visualization...');

for iFile = 2:3%length(imageFiles)
    
    currentFileName = imageFiles(iFile).name;
    
    % Load and Preprocess the image
    currentImage = imread(currentFileName);
    currentImage = mat2gray(currentImage);    % Convert to normalized double (0-1)
    currentImage = imadjust(currentImage);    % Adjust intensity

    [ySize, xSize] = size(currentImage);

    % Skip if the image is too small
    if ySize < subwindowSize(1) || xSize < subwindowSize(2)
        warning('Image %s is too small (%dx%d) for sub-window size %dx%d. Skipping.', ...
                currentFileName, ySize, xSize, subwindowSize(1), subwindowSize(2));
        continue; 
    end
    
    % Calculate dynamic step size and centers
    stepSize = subwindowSize .* stepFactor;
    %xCenters = stepSize(2):stepSize(2):(xSize - subwindowSize(2)/2);
    %yCenters = stepSize(1):stepSize(1):(ySize - subwindowSize(1)/2);
    half_window_y = subwindowSize(1)/2;
    half_window_x = subwindowSize(2)/2;

    xCenters = half_window_x : stepSize(2) : (xSize - half_window_x);
    yCenters = half_window_y : stepSize(1) : (ySize - half_window_y);
    
    if isempty(xCenters) || isempty(yCenters)
        warning('Image %s yielded no sub-windows. Skipping.', currentFileName);
        continue;
    end

    numWindowsX = length(xCenters);
    numWindowsY = length(yCenters);
    angles_xy = zeros(numWindowsY, numWindowsX);
    centerX_local = subwindowSize(2) / 2; % Local center for analyzeRotation
    centerY_local = subwindowSize(1) / 2; % Local center for analyzeRotation
    centerZ_local = 1;

    % --- First Pass & Second Pass Angle Calculation ---
    % Combined loop structure for both passes
    
    % 1. Run First Pass to fill initial angles_xy matrix
    for iX = 1:numWindowsX
        for iY = 1:numWindowsY
            xCenter_global = xCenters(iX);
            yCenter_global = yCenters(iY);
            
            xMin = round(xCenter_global - subwindowSize(2)/2 + 1);
            xMax = round(xCenter_global + subwindowSize(2)/2);
            yMin = round(yCenter_global - subwindowSize(1)/2 + 1);
            yMax = round(yCenter_global + subwindowSize(1)/2);
            subWindow = currentImage(yMin:yMax, xMin:xMax);
            
            %[minWidthAngle_xy, ~, ~, ~] = analyzeRotation(subWindow, centerX_local, centerY_local, centerZ_local, 'xy', angles, Rot_points, pixelpernm, false); % Set doVisualize to false here
            [minWidthAngle_xy, ~, ~] = manualFibrilOrientation(currentImage, subWindow, xMin, yMin, subwindowSize, currentFileName);
            angles_xy(iY, iX) = minWidthAngle_xy;
        end
    end
 
    %% 2. Run Second Pass (outlier check)
    %{
    padded_angles = NaN(numWindowsY + 2, numWindowsX + 2);
    padded_angles(2:end-1, 2:end-1) = angles_xy;
    
    for iX_pad = 2:numWindowsX + 1
        for iY_pad = 2:numWindowsY + 1
            xCenter_global = xCenters(iX_pad-1);
            yCenter_global = yCenters(iY_pad-1);
            
            neighbors = padded_angles(iY_pad-1:iY_pad+1, iX_pad-1:iX_pad+1);
            neighbors(2,2) = NaN;
            mean_neighbor = mean(neighbors(:),'omitnan');
            std_neighbor = std(neighbors(:),'omitnan');
            currentAngle = padded_angles(iY_pad, iX_pad);

            if abs(currentAngle - mean_neighbor) > std_neighbor
                % Recalculate range and get subWindow_pass2
                sweep_min = max(mean_neighbor - std_neighbor, min(angles));
                sweep_max = min(mean_neighbor + std_neighbor, max(angles));
                if (sweep_max - sweep_min) < 5, sweep_min = mean_neighbor - 2.5; sweep_max = mean_neighbor + 2.5; end
                sweep_min = max(sweep_min, min(angles)); sweep_max = min(sweep_max, max(angles));
                finerAngles = sweep_min:0.5:sweep_max;

                xMin = round(xCenter_global - subwindowSize(2)/2 + 1);
                xMax = round(xCenter_global + subwindowSize(2)/2);
                yMin = round(yCenter_global - subwindowSize(1)/2 + 1);
                yMax = round(yCenter_global + subwindowSize(1)/2);
                subWindow_pass2 = currentImage(yMin:yMax, xMin:xMax);

                % Re-calculate Rot_points for finerAngles
                half_size_finer = subwindowSize(1)/2;
                x1_finer = half_size_finer - (half_size_finer) * cosd(finerAngles'-90);
                y1_finer = half_size_finer - (half_size_finer) * sind(finerAngles'-90);
                x2_finer = half_size_finer + (half_size_finer) * cosd(finerAngles'-90);
                y2_finer = half_size_finer + (half_size_finer) * sind(finerAngles'-90);
                Rot_points_finer = [x1_finer, y1_finer, x2_finer, y2_finer];
                
                % Recalculate optimal angle
                [newAngle, ~, ~, ~] = analyzeRotation(subWindow_pass2, centerX_local, centerY_local, centerZ_local, 'xy', finerAngles, Rot_points_finer, pixelpernm, false);
                padded_angles(iY_pad, iX_pad) = newAngle;
            end
        end
    end
  
    angles_xy = padded_angles(2:end-1, 2:end-1);
  %}
    %% Step 5: Visualization of Fibril Directions
    
    % Define grid points for each sub-window center
    [xGrid, yGrid] = meshgrid(xCenters, yCenters);
    
    % Vector components (U, V) based on final angles
    U = cosd(angles_xy); % X-component
    V = sind(angles_xy); % Y-component

    % Create color index based on angle (mapping -60 to 60 deg across colormap)
    % Angle range is 120 degrees (-60 to 60). Normalize to [0, 1] then scale to [1, nColors].
    color_idx = round(((angles_xy + 60) / 120) * (size(customMap,1)-1)) + 1;
    color_idx = min(max(color_idx, 1), size(customMap,1)); % Ensure within bounds

    hFig = figure('Name', ['Fibril Map: ' currentFileName]);
    imshow(currentImage, []); % Display the image

    colormap(customMap); % Apply the custom colormap

    hold on
    for m = 1:numWindowsY
        for n = 1:numWindowsX   
            % Get color based on the angle
            color = customMap(color_idx(m, n), :);
            
            % Plot the arrow with quiver
            quiver(xGrid(m,n), yGrid(m,n),...
                   U(m,n), V(m,n),...
                   'LineWidth', 3, 'AutoScaleFactor', 15, 'Color', color, ...
                   'ShowArrowHead', 'off','Alignment','center'); 
        end
    end
    hold off;
    axis on; % Show axes for context
    title(['Fibril Orientation Map: ' currentFileName]);
    
    [~, name, ~] = fileparts(currentFileName);
    outputFigureName = [name '_FibrilMap.png'];
    saveas(hFig, outputFigureName);
    disp(['Saved visualization to: ', outputFigureName]);
    
    %% Step 6: Format and Store Results for CSV
    
    final_angles = angles_xy(:);
    x_centers = xGrid(:);
    y_centers = yGrid(:);
    
    % Calculate Unit Vector (w=0)
    u_vec = cosd(final_angles); 
    v_vec = sind(final_angles); 
    w_vec = zeros(size(final_angles)); 

    numMeasurements = length(final_angles);
    
    for k = 1:numMeasurements
        file_name = {currentFileName};
        coords_str = {sprintf('(%d, %d, 0)', round(x_centers(k)), round(y_centers(k)))};
        unit_vector = [u_vec(k), -v_vec(k), w_vec(k)];
        
        resultsData(resultRowIndex, :) = [file_name, coords_str, num2cell(unit_vector)];
        resultRowIndex = resultRowIndex + 1;
    end
end

%% Step 7: Write Output CSV
disp('Writing final results to CSV...');

T = cell2table(resultsData, ...
          'VariableNames', {'FileName', 'SubvolumeCenter_xyz', 'Vector_u', 'Vector_v', 'Vector_w'});

outputFileName = ['fibril_orientation_2d_results_manual.csv'];
writetable(T, outputFileName);

disp(['Successfully wrote results to: ', outputFileName]);