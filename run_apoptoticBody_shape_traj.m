% function:  run_apoptoticBody_shape
% This function finds the edges of a deforming apoptoic body for each frame
% in a cropped movie subset as the apoptoic body passes through and leaves 
% the stagnation point region.

% Updates 2021:
% - get_round_objects provides tr_filt_vis, list of potential AB particles
% after filter by number of frames in trajectory and visual inspection for
% in focus, clarity
% - Now can analyze the whole uncropped movie at once, no need to create
% cropped substacks manually in ImageJ first
% - replace confusing scope_mode input parameter for
% get_edge_rays_grad_subpix_ph_drdth_rave_filters.m (ONLY!)
% with 'maxgrad' or 'mingrad' to select option for defining the particle
% edge as draw rays out from the center
% - Subtract off background before set up mask and find edges 7/27/21
% - Automated quality control to determine 'good' ABs based on fitted
% ellipse area and orientation angle vs. time. August 2021

% Updates 2022:
% - 1/11/22: Set edge check figures to not be visible so that figures are not
% displayed. Final plotted frame figures are still saved. Code will run faster. 

% clear
clear AB_edge_data AB_edge_data_info
% close all
% clc

% Copyright (c) 2021, Joanna B. Dahl
% All rights reserved.
%
% Modified from scan_data_v4 from bending modulus code primarily written by
% Vivek Narsimhan for vesicle dumbbell instability 2016 Soft Matter paper
% and run_sphere_shape_history from SB3C 2019 abstract and poster and 
% run_apoptoticBody_shape from BBA GEN 2021 paper.
%
% Plot checks
% figure(3):    rough_AB after image binarization and filling holes
% figure(4):    original image
% figure(5):    pre-processed image (smoothing)
% figure(6):    rough cell outline with thresholding
% figure(7):    rough cell outline as a mask
% figure(8):    dilated rough cell outline mask
% figure(9):    in get_edge_rays_grad, intensity along ray Irad and grad
% figure(10):   in get_edge_rays_grad, original image with overlayed
%               detected edge
% figure(14):   in get_edge_rays_grad, radius check with ave
% figure(11):   area, perimeter subplot for current sphere
% figure(12):   circularity, deformation parameter for current sphere
% figure(13):   fitted ellipse axes, strain for current sphere

%% Set analysis parameters

% Image processing
smthR = 2; % Gaussian smoothing radius for image pre-processing, px
sens = 0.6; % binarization sensitivity for adaptive threshold, not used after 7/27/21 update
search_L = 50;  % px, half-length of square region centered at AB 
                % approximate location found from particle trajectory
                % analysis. Previously when manually cropping AB substacks,
                % 100 x 100 µm to 150 x 150 µm was a common size.
nbhd = 12;      % size of neighborhood of dilation in mask making

% gl_thrhld_mask = 85; % update 7/27/21. Lower threshold than trajectory code to capture closer to edge
gl_thrhld_mask = 180; % update 9/30/21. Increase threshold a bit to 150, then 180 to better capture black ABs for 8/23/21 experiments. White with black halo still ok for now
% Got value of 100 from playing around in ImageJ, eyeballing white ABs
% looking like big black solid circles and black ABs looking like a black
% halo (with binarization)
% JBD Data 2/UMB MGH ABs/2021.04.25_Gli36_IDH1_WT_mutant/Gli36_WT/xslot_Gli36_WT_750uLhr_z30um_mov_1
% See screenshots in UMB MGH ABs/analysis code for stiffness/ImageJ_thresholding

% Sensitivity comments:
 % sens = 0.615 works for many ABs from MGH, early 2021 experiments 
 % sens = 0.62 for faint ABs, vesicle-like small ABs
 % sens = 0.605 for big lumpy ABs, clumps of material >3um diameter,
 % out-of-focus faint ABs
 % sens = 0.645 for 2019.09.22 P1 visMC test/100x_visMC_600uLhr_SP_mov_3/
 %                      found here: /Volumes/JBD Data 2/UMB exosomes/
 % sens = 0.5 ok for alginate particles, but picks up small particles.

% Plot choices
iplot_bkgr = 1;     % Plot checks for background image update 7/27/21
iplot_mask = 1;     % For plot checks in mask creation section, logical
iplot_edge = 1;     % For plot checks in edge detection section, logical,
                    % obsolete after 7/28/21 update to dynamically set #
                    % frames to plot check (evenly spaced)
iplot_num  = 1;     % plot every X frames
figFont    = 14;       % font size of figure titles
ivis       = 1;     % Visibility of plot checks for edge detection, logical 1/11/22
                    % 0 = visible off, 1 = visible on

% QUALITY CONTROL parameters
% These limits on ellipse area, angle, and strain determine which ABs are
% save, which are labeled as good or bad, and which frames are thrown out
% due to shape uncertainty.
thrhld_A_scatter = 0.25;    % max tolerable scatter in std(A)/mean(A), for labeling suspicious, Default 0.25 works for ABs in initial testing
thrhld_strain_zero = 0.05;  % mean strain threshold below which strain is basically zero. Phi trend QC filter bypassed if phi close to zero. Default 0.05 from testing
thrhld_phi_trend = 0.6;     % threshold for % of phi angles of a certain sign to be considered a trend (positive or negative orientation angle). Defaul 0.6 or 60%. Frames not in trend are discarded.
thrhld_strain_range = 0.11; % threshold for range of ellipse strain, for labeling suspicious, Default 0.15
thrhld_ratio_stdMeanStr = 0.5;  % threshold for ratio of std to strain to assess scatter in strain, for labeling suspicious, Default 0.5
thrhld_postQC_numframes = 5;    % number of frames below which have low confidence in edge fitting (out of focus, unrealistic edges)

% Previous inputs in scan_data_v4, may use some of this later
% INPUTS:
%         COM = [x y] coordinates of the center of mass
%         folder = folder containing movie frames that we want to examine
%         startframe = first frame to analyze
%         endframe = last frame to analyze
%         output_filename = orig images + detected vesicle outline movie
%         filename as a string. NOTE: To save as avi movie file, no extension is
%         necessary. For saving as multipage tiff, make sure
%         output_filename has the format: '*.tif', thought this cabability
%         is not yet functional
% OPTION: FUNCTION CALL
% function [TBD] = run_sphere_edge_time_history(startframe, endframe, output_filename)

%% Load experimental image files from movies
% Load image files. Note the different syntax for PCs and Macs:
%   PC notation use backslash between directories.
%       datafolder = sprintf('Videos\\%s\\', folder);
%       srcFiles = dir( sprintf('%s\\*.tif', datafolder) );
%   Mac notation uses forward slash between directories.
%   Need to add '/' to end of datafolder name of PC filenames to work on Mac
% 
% % --- Update 10/28/21 Unnecessary if using run_case_ABshape_traj.m --- %
% % Specify datafolder
% disp('Select datafolder containing ORIGINAL movie frames:')
% datafolder = strcat(uigetdir('', 'Select datafolder containing movie frames:'),'/');
% % --- Update 10/28/21 Unnecessary if using run_case_ABshape_traj.m --- %

% Find the movie tiff files datafolder
 disp('Select datafolder containing ORIGINAL movie frames:')
 datafolder = strcat(uigetdir('', 'Select datafolder containing movie frames:'),'/');
% SHORT CUT: read in test movie in this folder
% datafolder = cd;
% datafolder in that datafolder
srcFiles = dir( sprintf('%s//*.tif', datafolder) );  % individiual tiffs
%srcFiles(1) = []; % had to do this for 750uLhr because my modified was and 1250uLhr
%also saved and I couldn't delete.
N = length(srcFiles);
% append separator
datafolder = strcat(datafolder,filesep);
% print data folder to screen
fprintf('\n%s\n',datafolder);

if N == 1       % multi-page tiff
    flag_multi = 1;     % turn on flag for multipage tiff source files
    mov_info = imfinfo(strcat(datafolder,srcFiles.name));
    numframes   = length(mov_info);   % number of frames
elseif N > 1    % individual tiff files (image sequence)
    numframes = N;      % number of frames
    flag_multi = 0;     % multipage tiff flag off
else
    disp('Problem reading in tiff files.')
    return
end

% for development and debugging, specify datafolder
% datafolder = '/Users/joannabechtel_preXP/Dropbox/Matlab codes/Image
% processing cell edge/'; % LAPTOP
% numframes = endframe - startframe + 1;    % if using start, end frames

disp(['No. frames:  ',num2str(numframes)])

%%  Read images into MATLAB 3D array - don't need after save bkgrdImg in traj subsection
% % update 7/27/21
% 
% for nn = 1:numframes
%     if flag_multi
%         % Get filename, and load image for multipage tiff movies
%         filename = strcat(datafolder, srcFiles.name); 
%         image = imread(filename,nn);    % Load image
%         Data_stack(:,:,nn) = image;
%     else
%         % Get filename, and load image for movies as individual tiff files
%         filename = strcat(datafolder, srcFiles(nn).name);
%         image = imread(filename);   % Load image
%     end
% end
% 
% [rows, columns, numSlices] = size(Data_stack);
% 
% if iplot_bkgr, figure, h = gcf; imshow(Data_stack(:,:,round(numSlices/2)),'DisplayRange',[]), figure(h), title('Checking movie import. This is a frame in the middle of the stack.'), end

%% Z projection - average image = background for this movie - - don't need after save bkgrdImg in traj subsection
% % run("Z Project...", "projection=[Average Intensity]");
% % update 7/27/21
% 
% bkgrdImg = zeros(rows, columns, class(Data_stack)); % Or whatever class you want.
% for col = 1 : columns
%     for row = 1 : rows
%         thisZVector = Data_stack(row, col, :);
%         meanValue = mean(thisZVector);
%         bkgrdImg(row, col) = meanValue;
%     end
% end
% disp('Background image found.')
% 
% if iplot_bkgr, figure, h = gcf; imshow(bkgrdImg,[]), figure(h), title('Background = ZProj Average'), end

%% Load trajectory information
% The trajectories are found from his comes from the function 
% run_threshold_get_AB_trajectories.m that find particles and their 
% trajectories on background-subtracted, smoothed, thresholded movies.

% Update 8/3/21: Automatically pick out the ABtracks mat file

% Find AB_shape data files in this datafolder
ABtrackFiles = dir( sprintf('%s/*ABtracks*.mat', datafolder) );  % look for shape mats
% Eliminate the one that has BACKUP in title if using older version of 
% run_threshold_get_AB_trajectories.m
for ii = 1:size(ABtrackFiles,1)
    if contains(ABtrackFiles(ii).name, 'BACKUP')
        ABtrackFiles(ii) = [];
    end
end

% Debugging print to screen
% disp('Check list of ABtrack mat files. Should be 1x1')
% ABtrackFiles_SIZE = size(ABtrackFiles) % check output. Should be 1x1

load(strcat(datafolder, ABtrackFiles.name))
fprintf('ABtracks mat-file: \t%s\n',ABtrackFiles.name);

% Previous version
% disp('Pick AB trajectory data file')
% [datafile, datafolder, filterindex] = uigetfile('*m', ...
%     'Pick AB trajectory data file', datafolder);
% 
% load(strcat(datafolder, datafile))
% fprintf('%s\n',datafile);

AB_IDs = unique(tr_filt_vis(:,4));
fprintf('No. AB trajectories:  \t%i\n',length(AB_IDs))

%% START AB EDGE DETECTION LOOP

for ind_ID = 1:length(AB_IDs) % loop over each AB trajectory
 %for ind_ID = 20:length(AB_IDs) % loop subset of AB trajectory data
    
    close all   % close all figures at start of iteration
    
    %% Trajectory info for this AB
    ID = AB_IDs(ind_ID);    % particle number

    % Pick out rows in tr_filt that correspond to AB with current ID
    AB_traj_info = tr_filt_vis( tr_filt_vis(:,4)==ID ,:); 
    x_tr = AB_traj_info(:,1);       % approx. AB x coordinate [px]
    y_tr = AB_traj_info(:,2);       % approx. AB y coordinate [px]
    frames_tr = AB_traj_info(:,3);  % trajectory frames
    
    % Status update printed to screen
    fprintf('AB%02d (%d out of %d):  Mask and edge detection over %i frames\n\n',...
        ID, ind_ID, length(AB_IDs), length(frames_tr));  
    
    % Update 7/28/21 dynamically set the number of frames to plot with edge
    if iplot_edge
        if length(frames_tr) <= 10
            iplot_num = 2;
        elseif (length(frames_tr) > 10) && (length(frames_tr) <= 20)
            iplot_num = 3;
        elseif (length(frames_tr) > 20) && (length(frames_tr) <= 30)
            iplot_num = 4;
        else
            iplot_num = round(length(frames_tr)/10);
        end
    end
    % end 7/28/21 update
    
%% MASKS Detect This Apoptotic Body's Masks for each frame in its trajectory
% Update 2/21: Using tr_filt_vis variable from get_round_objects, generate
% local mask for each frame in this AB trajectory
% Loop through each frame in this AB tracjectory's movie subset and find 
% the mask for each frame in the trajectory.

 % number of frames in this AB trajectory
numframes_ID = length(frames_tr);

for ii = 1:numframes_ID  % mask creation loop for AB ID
    
    % Display loop number
%     fprintf('Mask  detection for frame %i\n', frames_tr(ii));

    % Load current frame ii
    if flag_multi
        % Get filename, and load image for multipage tiff movies
        filename = strcat(datafolder, srcFiles.name); 
        image = imread(filename,frames_tr(ii));    % Load image
    else
        % Get filename, and load image for movies as individual tiff files
        filename = strcat(datafolder, srcFiles(frames_tr(ii)).name);
        image = imread(filename);   % Load image
    end
    
    % JBD plot check
    if iplot_mask
        if ~mod(ii,iplot_num)
            figOrig = figure(4); 
            imshow(image,[]); % true px intensity values
            title('Original Image'); 
    %         pos = get(gcf,'Position'); % set(gcf,'Position',[82 451 pos(3) pos(4)])
        end
    end
    % JBD plot check
    
    % JBD 7/27/21 update - subtract off background before find mask
    image_orig = image;
    image = imsubtract(image, bkgrdImg);
    
%   Preprocess image, currently Gaussian smoothing smthRxsmthR nbhd and sigma = 0.5
    image = preprocess_image(image, smthR);
    % JBD plot check: again to see effect of image pre-processing
    if iplot_mask
        if ~mod(ii,iplot_num)
            figure(5); imshow(image,[]); 
            title('Pre-processed image'); 
        end
    end
    % JBD plot check
    
    % Use AB trajectory information to create a mask for edge detection at
    % the approximate AB location.
    % 'Crop' to only search in region close around the AB.
    search_COM = round([x_tr(ii); y_tr(ii)]);   % crop search region center
    % rectangle vertex coordinates [xUL yUL; xUR yUR; xLR yLR; xLL yLL]
    search_vertices = [ search_COM(1)-search_L search_COM(2)-search_L; 
                        search_COM(1)-search_L search_COM(2)+search_L;
                        search_COM(1)+search_L search_COM(2)+search_L;
                        search_COM(1)+search_L search_COM(2)-search_L]; % global coord
    % Search crop region adjustments to fit into image pixel locations
    % (1) Set any negative pixel locations to 2, one row/column in from the
    % edge of the image. Setting search region up to the edge of the image
    % leads to detecting a boundary at this flush edge.
    search_vertices(search_vertices < 1) = 2;
    % (2) Set pixels larger than image dimensions to the max image dimension-1
    numpx = size(image, 2); numpy = size(image,1);
    search_vertices(search_vertices(:,1) > numpx,1) = numpx-1; % global coord
    search_vertices(search_vertices(:,2) > numpy,2) = numpy-1; % global coord
    
    % Pick out upper left corner of search region for easier local-global
    % pixel coordinate system transformations
    UL = search_vertices(1,:); % upper left x, y coord
    
    % plot check for crop search region
    if iplot_mask
        if ~mod(ii,iplot_num)
            figure(5); hold on; % plot on pre-processed image
            % reverse to (y,x) for plotting on image? No, not here!
            fill(search_vertices(:,1), search_vertices(:,2), 'c', 'FaceAlpha', 0.2)
        end
    end
    % Threshold and binarize the search region. Recall x,y indices switch
    % for images plotted with imshow
    search_region = image(search_vertices(1,2):search_vertices(2,2), ... % y, columns
                          search_vertices(2,1):search_vertices(3,1));    % x, rows
                      
    % CONVERSION: local (search_region) -> global (image) coordinates
    % x_global = x_local + (UL(1) - 1)
    % y_global = y_local + (UL(2) - 1)
   
    % Binarize image - turn to black and white with adaptive threshold
    % Result dependent on the Sensitivity parameter. 
    % Similar result obtained with BW = im2bw(search_region)
    % Likely need to fine-tune sensitivity for each AB movie.

%     % Sometimes edges of image can yield sharp gradients from artifacts in
%     % the image. So cut out the outer rows of pixels.
%     inset = 2;
%     BW = imbinarize(image(inset+1:end-inset,inset+1:end-inset),'adaptive','ForegroundPolarity','bright','Sensitivity',sens);
    
    % Note 7/9/21: The imbinarize command below seems to struggle with dark
    % apoptotic bodies and dark hydrogels. The thresholding in the
    % particle_tracking seems to work week for both light and dark
    % apoptotic bodies. Maybe use that image binarizing technique over this
    % imbinarize command. If the BW variable accurately shows boundaries,
    % then the rest of the code to fill holes and find edges should work
    % ok. The problem may be rooted in setting the 'ForegroundPolarity' to
    % 'bright'.
    % replaced by 7/27/21 update - START
%     BW = imbinarize(search_region,'adaptive','ForegroundPolarity',...
%             'bright','Sensitivity', sens); % replaced by 7/27/21 update
%     % Flip values of black and white for imfill function
%     BWc = imcomplement(BW);
%     % Fill holes. Likely just detected the outline of the particle/cell
%     rough_AB = imfill(BWc,'holes');
%     % iplot check for correct sensitivity setting, cropped FOV
%     if iplot_mask
%         if ~mod(ii,iplot_num)
%             figure(3)
%             subplot(2,1,1), imshow(search_region,'DisplayRange',[])
%             title('Search region')
%             subplot(2,1,2), imshow([BWc rough_AB]); 
%             title('Rough AB after binarization (pre, post imfill)'); 
%         end
%     end 
    % replaced by 7/27/21 update - END        
        
    % Update 7/27/21 Thresholding like run_threshold_get_AB_trajectories.m
    % imbinarize function needs a scaled threshold between 0 and 1
    % image must have the background subtracted
    if isa(image,'uint16')
        T = gl_thrhld_mask/(2^16-1);     % 16-bit images
    elseif isa(image, 'uint8')
        T = gl_thrhld_mask/(2^8-1);      % 8-bit images
    else
        disp('Unsupported image class. Try again')
    end
    BW = imbinarize(search_region, T);
    rough_AB = BW;
    if iplot_mask
        if ~mod(ii,iplot_num)
            figure(3)
            subplot(1,3,1), imshow(search_region,[])
            title('Search region')
            subplot(1,3,2), imshow(rough_AB); 
            title('Rough AB after binarization'); 
        end
    end 
    % End update 7/27/21
            
    
    % Find edges using bwboundaries. This function expects a binary image input:
    % Binary image where nonzero pixels (value 1 here) are the object and 0
    % pixels constitute the background.
    [B,L] = bwboundaries(rough_AB,'noholes'); % local coord
    
    bdy_COM = []; % clear and preallocate
    for k = 1:length(B)
        boundary = B{k};
        % Find COM of boundaries, convert to global coord system
        bdy_COM(k,:) = [mean(boundary(:,2))+(UL(1)-1) ...
            mean(boundary(:,1))+(UL(2)-1)]; % global coord, image coord.
        % need to convert to global coord to plot on image
%                 plot(boundary(:,2)+inset, boundary(:,1)+inset, 'w', 'LineWidth', 2)
        % bwboundaries now with MATLAB R2020b seems to store point coordinates in
        % image convention (1,1) origin in upper left
        % plot all boundaries on pre-processed image
%         if iplot_mask
%             if ~mod(ii,iplot_num)
%                 figure(6)
%                 imshow(image,[]); hold on
%                 title('B/W boundaries after thresholding')               
%                 plot(boundary(:,2)+(UL(1)-1), ... % swap y,x
%                      boundary(:,1)+(UL(2)-1), 'w.', 'LineWidth', 1) % swap
%             end
%         end 
    end

    if isempty(bdy_COM)
        mask = zeros(size(image));  % all black
        continue    % pass control to next iteration of this for loop
    end
    
    % UPDATE 2/21 - start
    % Assume the object with  COM close to COM of AB trajectory point
    % is the particle/cell of interest. Any other objects 
    % could be small particles, junk, wrong object, etc. Switch x, y 
    % columns so with plot function use syntax plot(X,Y)
    dx = bdy_COM(:,1) - x_tr(ii); % x_tr(ii) is the x coordinate of a particle trajectory in this frame
    dy = bdy_COM(:,2) - y_tr(ii); % y_tr(ii) is the x coordinate of a particle trajectory
    dist = sqrt(dx.^2 + dy.^2);  % distance from each boundary COM to approx AB COM
    [min_dist,min_ind] = min(dist);
    % Maybe also need a check that number of boundary points is >50 to
    % eliminate single pixel strays
    boundary = B{min_ind}; % local coord., y,x swapped coord.    
    % UPDATE 2/21 - end
    
    % Previously: Assume the object with the largest number of points in 
    % the B is the particle/cell of interest.
%     [max_size, max_index] = max(cellfun('size', B, 1)); % pick biggest boundary
%     boundary = B{max_index};
    
%     rough_bdry_pts = [boundary(:,2)+inset, boundary(:,1)+inset]; 
    % convert rough boundary points of AB to global coord.
    % rough_bdry_pts format [y x]. Opposite b/c turned into mask
    rough_bdry_pts = [boundary(:,2)+(UL(1)-1), ... % x coord. global
                      boundary(:,1)+(UL(2)-1)]; % y coord, no default inset
    % COM format [x_COM y_COM]
    rough_bdry_COM = [mean(rough_bdry_pts(:,1)) mean(rough_bdry_pts(:,2))];
    if iplot_mask
        if ~mod(ii,iplot_num)
            figure(4), hold on
            plot(rough_bdry_pts(:,1), rough_bdry_pts(:,2), 'w.')
            hold off
            figure(6)
            imshow(image,[]); hold on
%             plot(rough_bdry_COM(1), rough_bdry_COM(2),'or','MarkerFaceColor','r')
            plot(rough_bdry_pts(:,1), rough_bdry_pts(:,2), 'w.')
            hold off
            title('Fullsize image with AB rough edge')
            figure(3)
            subplot(1,3,3)
            imshow(search_region,[]), hold on
            plot(boundary(:,2), boundary(:,1), '.w'), hold off
            title('Search region with AB rough edge')
        end
    end
    
    % Dilate the rough AB boundary points to create the mask
    % Create a mask from boundary points
    mask = zeros(size(image));  % all black
    for k = 1:size(rough_bdry_pts,1)
        xm = rough_bdry_pts(k,1); 
        ym = rough_bdry_pts(k,2);
        mask(ym,xm) = 1;
    end
    mask = logical(mask);   % convert to logical
    
    if iplot_mask
        if ~mod(ii,iplot_num)
            figure(7), imshow(mask)
            title('Rough Apoptotic Body Outline, Pre-Mask')
        end
    end
    
    SE = strel('disk',nbhd);    % define structuring element, disk
    mask = imdilate(mask,SE);   % dilate rough cell outline
    if iplot_mask
        if ~mod(ii,iplot_num)
            figure(8), 
            subplot(1,2,1), imshow(mask);
            title(['Mask from dilated rough apoptoic body, Frame No. ',num2str(ii)],...
                'FontSize',figFont); 
            subplot(1,2,2), imshow(labeloverlay(image_orig, mask,'Transparency',0.9))
            figure(8), title('Mask Over Original Image')
        end
%         pause(0.5)
    end
    
    % Now we have a mask for the edge detection!
    mask_all{ii} = mask;
   
    
end  % mask creation loop

disp('Mask creation done!')

% close all

%% EDGES Detect This Apoptoic Body's Edges
% Ultizing the masks detected in the previous section, detect edges of
% apoptoic bodies by drawing rays from the center. The edge is defined to
% subpixel accuracy to be the location of the largest gradient.

counter = 0;

clear AB_edge_data AB_edge_data_info

% Preallocate
% AB_edge_data.x = NaN(360,numframes_ID);
% AB_edge_data.y = NaN(360,numframes_ID);
AB_edge_data.x = NaN(360,1);
AB_edge_data.y = NaN(360,1);
% Pplot = zeros(numframes_ID, 1);
% Aplot = zeros(numframes_ID, 1);
% aplot = zeros(numframes_ID, 1);
% bplot = zeros(numframes_ID, 1);
% phiplot = zeros(numframes_ID, 1);
% strainplot = zeros(numframes_ID, 1);

for ii = 1:numframes_ID  % edge detection loop for AB ID

% for ii = [1 2 4]  % edge detection loop for AB ID, skip some
% frames
    
    counter = counter + 1;
    
    % Display loop number
    fprintf('Edge detection for global frame %d, local frame %d (of %d)\n',...
        frames_tr(ii), ii, numframes_ID);
    
    % Load current frame ii
    if flag_multi
        % Get filename, and load image for multipage tiff movies
        filename = strcat(datafolder, srcFiles.name); 
        image = imread(filename,frames_tr(ii));    % Load image
    else
        % Get filename, and load image for movies as individual tiff files
        filename = strcat(datafolder, srcFiles(frames_tr(ii)).name);
        image = imread(filename);   % Load image
    end
    
    % JBD 7/27/21 update - subtract off background before find mask
    image_orig = image;
    image = imsubtract(image, bkgrdImg);
    
    % Preprocess image, currently Gaussian smoothing, currently radius set
    % by user, sigma = 0.5
    image = preprocess_image(image, smthR);
    
    % Select the appropriate initial mask from mask array created in the
    % previous section
    mask = cell2mat(mask_all(ii));
    
    if max(max(mask)) == 0      % if mask all black, then skip this frame
        continue
    end
    
    % Obtain mask center of mass and x,y extrema for this initial mask for 
    % get_edge_rays contour finding
    [center_of_mass, x_extrema, y_extrema] = get_mask_features(mask);
    
%     % see effect of image pre-processing and mask location
%     % NOTE: Moved down after get_edge_rays
%     if iplot_edge
% %         if ~mod(ii,iplot_num)
%             figOrig = figure(4); 
%             imshow(image,[]); % true px intensity values
%             title(['Original Image: AB ',num2str(ID),', Global Frame ',...
%                 num2str(frames_tr(ii)),', Local Frame ',num2str(ii)]); 
%             figure(5); ax1 = gca;
%             imshow(labeloverlay(image,mask,'Transparency',0.9),[]); 
%             title('Pre-processed image and mask overlay','FontSize',figFont); hold on
% %         end
%     end
%     
%     if iplot_edge
%         if ~mod(ii,iplot_num)
%             figure(8)   % mask plot
%             imshow(mask); hold on;
%             plot(center_of_mass(1), center_of_mass(2), 'oy')
%             fill([x_extrema(1) x_extrema(1) x_extrema(2) x_extrema(2)],...
%                  [y_extrema(1) y_extrema(2) y_extrema(2) y_extrema(1)],...
%                  'y', 'FaceAlpha', 0.2)
%             hold off
%             title(['Mask, COM, x,y extrema: Frame ',num2str(ii)],'FontSize',figFont)
%         end
%     end
    
    % Find (x,y) coordinates of edge from the maximum gradient along
    % rays drawn out from the center of the apoptoic body
    % UPDATE 7/9/21: edge_def_mode replaces scope_mode
    edge_def_mode = 'maxgrad';      % black edge to white halo, best for apoptotic bodies and alginate particles
%     edge_def_mode = 'mingrad';

%     scope_mode = 'brightfield';     % Don't use anymore. best for alginate particles, apoptotic bodes,edge = max(grad)
%     scope_mode = 'phase';         % Obsolete now. best for cells with junky interior, edge = min(grad)
    [x, y] = get_edge_rays_grad_subpix_ph_drdth_rave_filters(image, center_of_mass, x_extrema, y_extrema, mask, edge_def_mode);
    %  UPDATE 7/9/21 end
    
    % Update 7/28/21
    if sum(isnan(x)+isnan(y)) > 2     % if points are NaN, then skip this frame
        fprintf('(x,y) contour points are NaN. Frame %i skipped.\n',ii)
        continue
    end
    % Update 7/28/21
    
%     % Convert x, y coordinates of edge into (rho, theta) coordinates
%     % centered at the contour center of mass. Easier to calculate
%     % perimeter and area from these coordinates.
%     [rho, theta, x0, y0] = get_rho_theta(x0, y0);
%     % Find sphere/cell area, perimeter, circularity, deformation measure
%     % def = 0 is a circle, d > 0 is another shape
%     [P, A, c, def] = get_shape_metrics(rho, theta);

    % Find sphere/cell area, perimeter, circularity, deformation measure
    % from (x,y) coordinates
    [P, A, c, def] = get_shape_metrics_xy(x, y);
    
    % Edge plot check
    if iplot_edge
        if ~mod(ii,iplot_num) 
            % Visibility off to not display 1/11/22
            if ivis == 0  % if visibility set to off
                figOrig = figure('visible','off');
            else
                figOrig = figure(4);
            end
                
            imshow(image_orig,[]); % true px intensity values
            title(['Original Image: AB ',num2str(ID),', Global Frame ',...
                num2str(frames_tr(ii)),', Local Frame ',num2str(ii),...
                ' of ',num2str(numframes_ID)]); 
            %set(gcf,'Position', [35 215 930 871])
            set(gcf,'Position', [35 215 530 471])
            
            % Visibility off to not display 1/11/22
            if ivis == 0  % if visibility set to off
                hfig5 = figure('visible','off');
            else
                hfig5 = figure(5);
            end
            
            ax1 = gca;
            imshow(image_orig,[]); hold on % true px intensity values
    %         imshow(labeloverlay(image,mask,'Transparency',0.9),[]); hold on
            plot(x,y,'.y')
            title(['Contour points: AB ',num2str(ID),', Global Frame ',num2str(frames_tr(ii)),...
                ', Local Frame ',num2str(ii),' of ',num2str(numframes_ID)]); 
            %set(gcf,'Position', [975 215 930 871])
            set(gcf,'Position', [575 215 530 471])
            
            if ivis % if figures visible, then pause
                pause(0.1)
            end
            
        else
            ax1 = [];   % send empty axis handle to fit_ellipse
        end

    end
    
    % Fit ellipse to contour
    % Find best fit ellipse and plot fitted ellipse on image in
    % figure(5) along with the detected contour points
    ellipse_t = fit_ellipse(x, y, ax1);
    
    % And calculate the the engineering strain 
    % strain = (a-b)/(a+b) where a is the long axis and b is the short
    % axis.
    if isempty(ellipse_t) %if ellipse fitting program issue (just one repeated contour point for instance)
        warning('No ellipse structure from fit_ellipse.m. Frame skipped.')
        continue
    elseif isempty(ellipse_t.a) % no ellipse fitted
        % something
        warning('No ellipse fitted within fit_ellipse.m. Frame skipped.')
        continue
    else
        a = ellipse_t.long_axis; b = ellipse_t.short_axis;  
        strain = (a - b)/(b + a);
        length_vert = sqrt(diff(ellipse_t.ver_line(:,1))^2 + diff(ellipse_t.ver_line(:,2))^2);
        length_horz = sqrt(diff(ellipse_t.horz_line(:,1))^2 + diff(ellipse_t.horz_line(:,2))^2);
    end
    
%     Pplot(counter) = P;
%     Aplot(counter) = A;
%     aplot(counter) = a;
%     bplot(counter) = b;
%     phiplot(counter) = ellipse_t.phi*180/pi;
%     strainplot(counter) = strain;
    
    % Save contour points, metrics, and other meta data in a structure array
    AB_edge_data.x(1:length(x),counter) = x;
    AB_edge_data.y(1:length(y),counter) = y;
    AB_edge_data.P(counter) = P;
    AB_edge_data.A(counter) = A;
    AB_edge_data.c(counter) = c;
    AB_edge_data.def(counter) = def;
    AB_edge_data.COM(:,counter) = [mean(x); mean(y)];
    AB_edge_data.a(counter) = a;
    AB_edge_data.b(counter) = b;
    AB_edge_data.phi(counter) = ellipse_t.phi;
    AB_edge_data.strain(counter) = strain;
    AB_edge_data.length_vert(counter) = length_vert;
    AB_edge_data.length_horz(counter) = length_horz;
    AB_edge_data.ellipse(counter) = ellipse_t;
    
end  % edge detection loop

disp('Edge detection done!')
% Output status to screen
fprintf('\nEdge detection, ellipse fitting complete for AB%02d (%d out of %d).\n\n',...
    ID,ind_ID,length(AB_IDs))

preQC_AB_edge_data = AB_edge_data;  % save a copy before start throwing out frames

%% Skip AB save for poor shape result, Remove frames with zero area
% Skip save if 
% - No ellipses for the AB
% Remove frames with zero area (no ellipse found

% No ellipses stored: skip plotting, saving
if ~isfield(AB_edge_data,'a')   
    warning('SKIP: No ellipses fitted to this AB. AB save skipped.')
    continue
end

% Only 1-2 frames with an AB ellipse shape, skip plotting, saving
% NOTE: Does not catch most events where plots only have 1-2 points.
% test output until figure out why not catching cases with less than 3
% plotted frames
length_a_AB_edge_data = length(AB_edge_data.a);
if length(AB_edge_data.a)<3   
    warning('SKIP: Fewer than 3 frames in shape analysis this AB. AB save skipped.')
    continue
end

% Delete frames with ellipse area of zero
% If get here, means at least one ellipse was detected - does it?
% Logical for frames wtih area == 0
log_zeroA = (AB_edge_data.A == 0);
if sum(log_zeroA) > 0
    fprintf('Eliminating %i frames with ellipse are of zero (no ellipse detected)\n',sum(log_zeroA));
    AB_edge_data.x(:,log_zeroA) = [];
    AB_edge_data.y(:,log_zeroA) = [];
    AB_edge_data.P(log_zeroA) = [];
    AB_edge_data.A(log_zeroA) = [];
    AB_edge_data.c(log_zeroA) = [];
    AB_edge_data.def(log_zeroA) = [];
    AB_edge_data.COM(:,log_zeroA) = [];
    AB_edge_data.a(log_zeroA) = [];
    AB_edge_data.b(log_zeroA) = [];
    AB_edge_data.phi(log_zeroA) = [];
    AB_edge_data.strain(log_zeroA) = [];
    AB_edge_data.length_vert(log_zeroA) = [];
    AB_edge_data.length_horz(log_zeroA) = [];
    AB_edge_data.ellipse(log_zeroA) = [];
end

    
%% Plot shape measurements - pre QC

% Visibility off to not display 1/11/22
if ivis  % if visibility set to off, don't bother with these plots

% PERIMETER
hfig(1) = figure(11);
subplot(2,1,1)
plot(AB_edge_data.P,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Perimeter [px]','FontSize',figFont)
title(['Sphere number ',num2str(ii)],'FontSize',figFont)
grid on
subplot(2,1,2)
plot(AB_edge_data.A,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Area [px^2]','FontSize',figFont)
set(gcf,'Position',[450  545   560   800])
grid on

% figure(12)
% subplot(2,1,1)
% plot(AB_edge_data(:).c,'.k')
% xlabel('Local frame number','FontSize',figFont) 
% ylabel('Circularity [-]','FontSize',figFont)
% title(['Sphere number ',num2str(ii)],'FontSize',figFont)
% subplot(2,1,2)
% plot(AB_edge_data(:).def,'.k')
% xlabel('Local frame number','FontSize',figFont) 
% ylabel('Deformation parameter d=1-c [-]','FontSize',figFont)

% ELLIPSE AXES, TILT ANGLE, STRAIN
hfig(2) = figure(13);
subplot(3,1,1)
plot(AB_edge_data.a,'or','MarkerFaceColor','r'); hold on
plot(AB_edge_data.b,'ob','MarkerFaceColor','b'); hold off
xlabel('Local frame number','FontSize',figFont) 
ylabel('Fitted Ellipse Axes [px]','FontSize',figFont)
legend('Long','Short')
title(['Sphere number ',num2str(ii)],'FontSize',figFont)
grid on
subplot(3,1,2)
plot(AB_edge_data.phi*180/pi,'ob','MarkerFaceColor','b'); % degrees
xlabel('Local frame number','FontSize',figFont) 
ylabel('Fitted Ellipse Orientation [degrees]','FontSize',figFont)
grid on
subplot(3,1,3)
plot(AB_edge_data.strain,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Strain \epsilon = (a-b)/(a+b) [-]','FontSize',figFont)
set(gcf,'Position',[1000  545   560   800])
grid on

end % end ivis don't plot

%% QUALITY CONTROL: Assess strain scatter, range of ellipse area and angle
% Skip AB save if scatter in area and angle is too lar

% Area population statistics
A_mean = mean(AB_edge_data.A);
A_std  = std(AB_edge_data.A);
A_max  = max(AB_edge_data.A);
A_min  = min(AB_edge_data.A);

% Linear regression of area with local frame number
[slope_A, int_A, s2_A, ~, ~, ~, ~,...
    r2_A, res_A, res_std_A] = linreg_lsq(1:length(AB_edge_data.A), AB_edge_data.A);  

% Get rid of area outliers here. Otherwise can delete a good AB if several
% bad frames out of tens of frames.
% Throw out frames with standardized residual absolute value > 2
outlier_A = abs(res_std_A)>2;
if sum(outlier_A) > 0
    fprintf('Eliminating %i frames with outlier areas (linear regression), initial round\n',sum(outlier_A));
    AB_edge_data.x(:,outlier_A) = [];
    AB_edge_data.y(:,outlier_A) = [];
    AB_edge_data.P(outlier_A) = [];
    AB_edge_data.A(outlier_A) = [];
    AB_edge_data.c(outlier_A) = [];
    AB_edge_data.def(outlier_A) = [];
    AB_edge_data.COM(:,outlier_A) = [];
    AB_edge_data.a(outlier_A) = [];
    AB_edge_data.b(outlier_A) = [];
    AB_edge_data.phi(outlier_A) = [];
    AB_edge_data.strain(outlier_A) = [];
    AB_edge_data.length_vert(outlier_A) = [];
    AB_edge_data.length_horz(outlier_A) = [];
    AB_edge_data.ellipse(outlier_A) = [];
    % Update linear regression
    [slope_A, int_A, s2_A, se_slope, se_int, conf95_slope, conf95_int,...
        r2_A, res_A, res_std_A] = linreg_lsq(1:length(AB_edge_data.A), AB_edge_data.A);  
end


%% QUALITY CONTROL: Ellipse angle
% -	Throw out frames for which angle sign (+/-) is opposite the main trend, 
% except if oscillating around zero or strain is small and possibly if frames 
% are close to evenly split +/- angles.
% - Main trend defined as >=60% of frames with that ellipse angle sign
% -	Consider zero angle to be ± 15˚ (or 10˚?)
% -	If AB strain < 0.05, don’t worry about angle flipping between positive 
% and negative values. The shape is very close to spherical, and these 
% small strains are hard to resolve.

% Angle population statistics, convert to degrees
phi_mean = mean(AB_edge_data.phi)*180/pi;
phi_std  = std(AB_edge_data.phi)*180/pi;
phi_max  = max(AB_edge_data.phi)*180/pi;
phi_min  = min(AB_edge_data.phi)*180/pi;

% Strain mean. Small strains mean angle can easily vary between +/-, even
% for clear edges
strain_mean = mean(AB_edge_data.strain);

% If angles are close to zero, likely scattered somewhat evenly around
% zero. Bypass further angle-based QC for this AB
% Definition of close to zero:
%   - Mean is close to zero: +/- 15˚
%   - Scattered around zero with tight-ish range: Assess with std
%   and/or range
phi_range = abs(phi_max - phi_min);  % degrees
if (phi_range < 45) && (abs(phi_mean) <=10)
    
    fprintf('Ellipse angle close to zero. Bypass angle-based QC.\n')
    phi_nonzero = 0;
    phi_suspicious = 0; % not suspicous if phi close to zero
    
else
    
    phi_nonzero = 1;
    % Determine main ellipse angle trend, if strains are large enough and phi
    % is non-zero. Initial testing used strain_mean > 0.05
    if (strain_mean > thrhld_strain_zero) && phi_nonzero

        % Determine number of positive and negative ellipse orientation angles in
        % movie frames, also as a percentage of total frames
        phi_pos = sum(AB_edge_data.phi >= 0);
        phi_neg = sum(AB_edge_data.phi < 0);
        per_pos = phi_pos/length(AB_edge_data.phi);
        per_neg = phi_neg/length(AB_edge_data.phi);

        % Trend = more than thrhld_phi_trend (60% default) of frames have 
        % ellipses angled in a given direction
        if (per_pos > per_neg) && (per_pos >= thrhld_phi_trend)
            fprintf('Mostly positive ellipse orientation angles. Discarding frames with phi < 0.\n');
            fprintf('Eliminating %i frames with phi < 0\n', phi_neg);
            % Positive ellipse angle is the dominant trend. Throw out frames
            % with negative angles.
            neg_phi_ind = AB_edge_data.phi < 0;     % indices of negative angle frames
            AB_edge_data.x(:,neg_phi_ind) = [];
            AB_edge_data.y(:,neg_phi_ind) = [];
            AB_edge_data.P(neg_phi_ind) = [];
            AB_edge_data.A(neg_phi_ind) = [];
            AB_edge_data.c(neg_phi_ind) = [];
            AB_edge_data.def(neg_phi_ind) = [];
            AB_edge_data.COM(:,neg_phi_ind) = [];
            AB_edge_data.a(neg_phi_ind) = [];
            AB_edge_data.b(neg_phi_ind) = [];
            AB_edge_data.phi(neg_phi_ind) = [];
            AB_edge_data.strain(neg_phi_ind) = [];
            AB_edge_data.length_vert(neg_phi_ind) = [];
            AB_edge_data.length_horz(neg_phi_ind) = [];
            AB_edge_data.ellipse(neg_phi_ind) = [];

            phi_suspicious = 0;  % trend is not suspicous
    
        elseif (per_neg > per_pos) && (per_neg >= thrhld_phi_trend)
            fprintf('Mostly negative ellipse orientation angles. Discarding frames with phi > 0.\n');
            fprintf('Eliminating %i frames with phi > 0\n', phi_pos);
            % Negative ellipse angle is the dominant trend. Throw out frames
            % with positive angles.
            pos_phi_ind = AB_edge_data.phi >= 0;     % indices of negative angle frames
            AB_edge_data.x(:,pos_phi_ind) = [];
            AB_edge_data.y(:,pos_phi_ind) = [];
            AB_edge_data.P(pos_phi_ind) = [];
            AB_edge_data.A(pos_phi_ind) = [];
            AB_edge_data.c(pos_phi_ind) = [];
            AB_edge_data.def(pos_phi_ind) = [];
            AB_edge_data.COM(:,pos_phi_ind) = [];
            AB_edge_data.a(pos_phi_ind) = [];
            AB_edge_data.b(pos_phi_ind) = [];
            AB_edge_data.phi(pos_phi_ind) = [];
            AB_edge_data.strain(pos_phi_ind) = [];
            AB_edge_data.length_vert(pos_phi_ind) = [];
            AB_edge_data.length_horz(pos_phi_ind) = [];
            AB_edge_data.ellipse(pos_phi_ind) = [];

            phi_suspicious = 0;  % trend is not suspicous

        else

            % Suspect ABs:
            % 1. Number of pos, negative angle frames are the same or no dominant
            % trend. 
            % 2. Phi range significant, >45 degrees so classified as nonzero
            warning('SUSPICIOUS phi angle: sizable range and even pos, neg angles.')
            phi_suspicious = 1;

        end

    else    % Average phi close to zero, but has wide range
        % No frames thrown out, but labeled as suspicious
        warning('SUSPICIOUS phi angle: average angle close to zero, but sizable range.')
        phi_suspicious = 1;

    end
    
end


%% QUALITY CONTROL: Area outliers

% Linear regression of area with local frame number
% [linCoeff_A, S] = polyfit(1:length(AB_edge_data.A), AB_edge_data.A, 1);
[slope_A, int_A, s2_A, se_slope, se_int, conf95_slope, conf95_int,...
    r2_A, res_A, res_std_A] = linreg_lsq(1:length(AB_edge_data.A), AB_edge_data.A);  

% Throw out frames with standardized residual absolute value > 2
outlier_A = abs(res_std_A)>2;
if sum(outlier_A) > 0
    fprintf('Eliminating %i frames with outlier areas (linear regression), after phi trend QC\n',sum(outlier_A));
    AB_edge_data.x(:,outlier_A) = [];
    AB_edge_data.y(:,outlier_A) = [];
    AB_edge_data.P(outlier_A) = [];
    AB_edge_data.A(outlier_A) = [];
    AB_edge_data.c(outlier_A) = [];
    AB_edge_data.def(outlier_A) = [];
    AB_edge_data.COM(:,outlier_A) = [];
    AB_edge_data.a(outlier_A) = [];
    AB_edge_data.b(outlier_A) = [];
    AB_edge_data.phi(outlier_A) = [];
    AB_edge_data.strain(outlier_A) = [];
    AB_edge_data.length_vert(outlier_A) = [];
    AB_edge_data.length_horz(outlier_A) = [];
    AB_edge_data.ellipse(outlier_A) = [];
    % Update linear regression
    [slope_A, int_A, s2_A, se_slope, se_int, conf95_slope, conf95_int,...
        r2_A, res_A, res_std_A] = linreg_lsq(1:length(AB_edge_data.A), AB_edge_data.A);  
end

% Area scatter considered too big if std >25% of the magnitude of the mean area
% AND standard deviation of linear regression >25% of mean area
% This means that steady linear trends in ellipse area due to vertical
% drift are allowed, even if the range is large.
% Testing reveals thrhld_A_scatter = 0.25 is pretty good for ABs
if (A_std/A_mean > thrhld_A_scatter) && (sqrt(s2_A)/A_mean > thrhld_A_scatter)
    % why checking before throw out zero area frames (A_std) and after sqrt(s2_A)?
    % I forget, but keeping for now.
    warning('SUSPICIOUS area: Ellipse area is suspiciously scattered over a wide range.')
    A_suspicious = 1;
else
    % Area is sufficiently tight
    A_suspicious = 0;
end

%% Plot check - mid QC

% Visibility off to not display 1/11/22
if ivis  % if visibility set to off, don't bother with these plots

% perimeter
hfig(1) = figure(11);
subplot(2,1,1)
plot(AB_edge_data.P,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Perimeter [px]','FontSize',figFont)
title(['Sphere number ',num2str(ii)],'FontSize',figFont)
grid on
subplot(2,1,2)
plot(AB_edge_data.A,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Area [px^2]','FontSize',figFont)
set(gcf,'Position',[450  545   560   800])
grid on

% figure(12)
% subplot(2,1,1)
% plot(AB_edge_data(:).c,'.k')
% xlabel('Local frame number','FontSize',figFont) 
% ylabel('Circularity [-]','FontSize',figFont)
% title(['Sphere number ',num2str(ii)],'FontSize',figFont)
% subplot(2,1,2)
% plot(AB_edge_data(:).def,'.k')
% xlabel('Local frame number','FontSize',figFont) 
% ylabel('Deformation parameter d=1-c [-]','FontSize',figFont)

hfig(2) = figure(13);
subplot(3,1,1)
plot(AB_edge_data.a,'or','MarkerFaceColor','r'); hold on
plot(AB_edge_data.b,'ob','MarkerFaceColor','b'); hold off
xlabel('Local frame number','FontSize',figFont) 
ylabel('Fitted Ellipse Axes [px]','FontSize',figFont)
legend('Long','Short')
title(['Sphere number ',num2str(ii)],'FontSize',figFont)
grid on
subplot(3,1,2)
plot(AB_edge_data.phi*180/pi,'ob','MarkerFaceColor','b'); % degrees
xlabel('Local frame number','FontSize',figFont) 
ylabel('Fitted Ellipse Orientation [degrees]','FontSize',figFont)
grid on
subplot(3,1,3)
plot(AB_edge_data.strain,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Strain \epsilon = (a-b)/(a+b) [-]','FontSize',figFont)
set(gcf,'Position',[1000  545   560   800])
grid on

end % end ivis don't plot

%% QUALITY CONTROL: Strain
% Strain: Ideally scattered around a mean value, even with vertical drift. 
% Label AB ‘bad’ if strain increases/decreases too much or has a big range 
% (>0.15 for now).

% Linear regression of strain with local frame number
[slope_strain, int_strain, s2_strain, se_slope_strain, se_int_strain,...
    conf95_slope_strain, conf95_int_strain,r2_strain, res_strain, res_std_strain] ...
    = linreg_lsq(1:length(AB_edge_data.strain), AB_edge_data.strain);
% Outliers - `
outlier_strain = abs(res_std_strain)>2;
if sum(outlier_strain) > 0
    fprintf('Eliminating %i frames with outlier strain (linear regression), after area outlier, phi trend QC\n',sum(outlier_strain));
    AB_edge_data.x(:,outlier_strain) = [];
    AB_edge_data.y(:,outlier_strain) = [];
    AB_edge_data.P(outlier_strain) = [];
    AB_edge_data.A(outlier_strain) = [];
    AB_edge_data.c(outlier_strain) = [];
    AB_edge_data.def(outlier_strain) = [];
    AB_edge_data.COM(:,outlier_strain) = [];
    AB_edge_data.a(outlier_strain) = [];
    AB_edge_data.b(outlier_strain) = [];
    AB_edge_data.phi(outlier_strain) = [];
    AB_edge_data.strain(outlier_strain) = [];
    AB_edge_data.length_vert(outlier_strain) = [];
    AB_edge_data.length_horz(outlier_strain) = [];
    AB_edge_data.ellipse(outlier_strain) = [];
    % Update linear regression
    [slope_strain, int_strain, s2_strain, se_slope_strain, se_int_strain,...
        conf95_slope_strain, conf95_int_strain,r2_strain, res_strain, res_std_strain] ...
        = linreg_lsq(1:length(AB_edge_data.strain), AB_edge_data.strain);
    outlier2_strain = abs(res_std_strain)>2;
end

% strain population statistics
strain_mean = mean(AB_edge_data.strain);
strain_std  = std(AB_edge_data.strain);
strain_max  = max(AB_edge_data.strain);
strain_min  = min(AB_edge_data.strain);
strain_range = abs(strain_max - strain_min);

% Suspicous strain measurements have large std compared to mean and have a
% wide range. Defaults thrhld_strain_range = 0.15, thrhld_ratio_stdMeanStr
% = 0.5
if (strain_range > thrhld_strain_range) || (strain_std/strain_mean > thrhld_ratio_stdMeanStr)
    strain_suspicious = 1;
    warning('SUSPICIOUS strain: sizable range or relatively large std.')
else
    strain_suspicious = 0;
end

%% QUALITY CONTROL: Number of frames

if length(AB_edge_data.strain) < thrhld_postQC_numframes % if less than minimum number of frames
    numFr_suspicious = 1;
    warning('SUSPICIOUS number of frames: Small number of frames.')
else
    numFr_suspicious = 0;
end

%% Plot check - post QC

% PERIMETER
% Visibility off to not display 1/11/22
if ivis == 0  % if visibility set to off
    hfig(1) = figure('visible','off');
else
    hfig(1) = figure(11);
end

subplot(2,1,1)
plot(AB_edge_data.P,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Perimeter [px]','FontSize',figFont)
title(['AB ',num2str(ID)],'FontSize',figFont);
grid on
subplot(2,1,2)
plot(AB_edge_data.A,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Area [px^2]','FontSize',figFont)
set(gcf,'Position',[450  545   560   800])
grid on

% figure(12)
% subplot(2,1,1)
% plot(AB_edge_data(:).c,'.k')
% xlabel('Local frame number','FontSize',figFont) 
% ylabel('Circularity [-]','FontSize',figFont)
% title(['Sphere number ',num2str(ii)],'FontSize',figFont)
% subplot(2,1,2)
% plot(AB_edge_data(:).def,'.k')
% xlabel('Local frame number','FontSize',figFont) 
% ylabel('Deformation parameter d=1-c [-]','FontSize',figFont)

% ELLIPSE AXES, TILT ANGLE, STRAIN
% Visibility off to not display 1/11/22
if ivis == 0  % if visibility set to off
    hfig(2) = figure('visible','off');
else
    hfig(2) = figure(13);
end

subplot(3,1,1)
plot(AB_edge_data.a,'or','MarkerFaceColor','r'); hold on
plot(AB_edge_data.b,'ob','MarkerFaceColor','b'); hold off
xlabel('Local frame number','FontSize',figFont) 
ylabel('Fitted Ellipse Axes [px]','FontSize',figFont)
legend('Long','Short')
title(['AB ',num2str(ID)],'FontSize',figFont);
grid on
subplot(3,1,2)
plot(AB_edge_data.phi*180/pi,'ob','MarkerFaceColor','b'); % degrees
xlabel('Local frame number','FontSize',figFont) 
ylabel('Fitted Ellipse Orientation [degrees]','FontSize',figFont)
grid on
subplot(3,1,3)
plot(AB_edge_data.strain,'ok','MarkerFaceColor','k')
xlabel('Local frame number','FontSize',figFont) 
ylabel('Strain \epsilon = (a-b)/(a+b) [-]','FontSize',figFont)
set(gcf,'Position',[1000  545   560   800])
grid on

%% Edge data MAT file name
% Contains judgement of edge quality and strain certainty
% Angle phi is usually suspicious. Give that one less weight
sum_suspicious = A_suspicious + phi_suspicious + strain_suspicious + numFr_suspicious;

fprintf('Suspicion metrics\n');
fprintf('Area: %i\t\t Phi: %i\t\t Strain: %i\t Num Fr: %i\n', A_suspicious,...
    phi_suspicious, strain_suspicious, numFr_suspicious);

if sum_suspicious >= 3
    % If 3+ flags, then label as 'bad'
    qual = 'bad';
elseif (sum_suspicious == 2) && (phi_suspicious == 0)
    % usually bad if 2 other than phi are suspicious
    qual = 'bad';
elseif (sum_suspicious == 2) && (phi_suspicious == 1)
    % phi is usually suspicious and should not tip the assessment to bad,
    % just ok
    qual = 'ok';
elseif sum_suspicious == 1
    % just one, then label 'ok'
    qual = 'ok';
else 
    % nothing suspcious then label 'good'
    qual = 'good';
end


%% OBSOLETE Delete spurious frames - final manual check
% Manual and time consuming. Replaced by auto quality control above

% flag_del = input('Delete any frames?  1 = yes, 0 = no      ');
flag_del = 0;

if flag_del
    frame_num = input('Frame number(s) to delete (as a vector):   ');

    AB_edge_data.x(:,frame_num) = [];
    AB_edge_data.y(:,frame_num) = [];
    AB_edge_data.P(frame_num) = [];
    AB_edge_data.A(frame_num) = [];
    AB_edge_data.c(frame_num) = [];
    AB_edge_data.def(frame_num) = [];
    AB_edge_data.COM(:,frame_num) = [];
    AB_edge_data.a(frame_num) = [];
    AB_edge_data.b(frame_num) = [];
    AB_edge_data.phi(frame_num) = [];
    AB_edge_data.strain(frame_num) = [];
    AB_edge_data.length_vert(frame_num) = [];
    AB_edge_data.length_horz(frame_num) = [];
    AB_edge_data.ellipse(frame_num) = [];
    
    % replot perimeter
    hfig(1) = figure(11);
    subplot(2,1,1)
    plot(AB_edge_data.P,'ok','MarkerFaceColor','k')
    xlabel('Local frame number','FontSize',figFont) 
    ylabel('Perimeter [px]','FontSize',figFont)
    title(['Sphere number ',num2str(ii)],'FontSize',figFont)
    grid on
    subplot(2,1,2)
    plot(AB_edge_data.A,'ok','MarkerFaceColor','k')
    xlabel('Local frame number','FontSize',figFont) 
    ylabel('Area [px^2]','FontSize',figFont)
    set(gcf,'Position',[450  545   560   800])
    grid on

    % figure(12)
    % subplot(2,1,1)
    % plot(AB_edge_data.c,'.k')
    % xlabel('Local frame number','FontSize',figFont) 
    % ylabel('Circularity [-]','FontSize',figFont)
    % title(['Sphere number ',num2str(ii)],'FontSize',figFont)
    % subplot(2,1,2)
    % plot(AB_edge_data.def,'.k')
    % xlabel('Local frame number','FontSize',figFont) 
    % ylabel('Deformation parameter d=1-c [-]','FontSize',figFont)

    hfig(2) = figure(13);
    subplot(3,1,1)
    plot(AB_edge_data.a,'or','MarkerFaceColor','r'); hold on
    plot(AB_edge_data.b,'ob','MarkerFaceColor','b'); hold off
    xlabel('Local frame number','FontSize',figFont) 
    ylabel('Fitted Ellipse Axes [px]','FontSize',figFont)
    legend('Long','Short')
    title(['Sphere number ',num2str(ii)],'FontSize',figFont)
    grid on
    subplot(3,1,2)
    plot(AB_edge_data.phi*180/pi,'ob','MarkerFaceColor','b'); % degrees
    xlabel('Local frame number','FontSize',figFont) 
    ylabel('Fitted Ellipse Orientation [degrees]','FontSize',figFont)
    grid on
    subplot(3,1,3)
    plot(AB_edge_data.strain,'ok','MarkerFaceColor','k')
    xlabel('Local frame number','FontSize',figFont) 
    ylabel('Strain \epsilon = (a-b)/(a+b) [-]','FontSize',figFont)
    set(gcf,'Position',[1000  545   560   800])
    grid on
end

%% Save data

length_a_AB_edge_data = length(AB_edge_data.a);
if length(AB_edge_data.a)<3   
    warning('SKIP: Post QC, fewer than 3 frames in shape analysis this AB. AB save skipped.')
    continue
end

% flag_save = input('Save the apoptoic body edge data structure? 1 = yes, 0 = no      ');
flag_save = 1;

if flag_save
    % New structure for metadata
    AB_edge_data_info.datafolder = datafolder;  % movie file location  
    AB_edge_data_info.movie_name = srcFiles.name; % movie file name
    AB_edge_data_info.smthR = smthR; % Gaussian smoothing radius
    AB_edge_data_info.sens = sens;  % auto threshold sensitivity
    AB_edge_data_info.gl_thrhld_mask = gl_thrhld_mask; % binarization thresholding
    AB_edge_data_info.search_L = search_L;  % px, half-length of square region centered at AB          
    AB_edge_data_info.nbhd = nbhd;  % px, size of neighborhood of dilation in mask making
    AB_edge_data_info.ID_tr = ID;    % AB ID from tr_filt_vis
    AB_edge_data_info.x_tr = x_tr;  % approx. AB x coordinate [px] from tr_filt_vis
    AB_edge_data_info.y_tr = y_tr;  % approx. AB y coordinate [px] from tr_filt_vis
    AB_edge_data_info.frames_tr = frames_tr;  % trajectory frames from tr_filt_vis
    AB_edge_data_info.thrhld_A_scatter = thrhld_A_scatter;  % QC threshold std(A)/mean(A), to label suspicious
    AB_edge_data_info.thrhld_strain_zero = thrhld_strain_zero; % QC threshold below which strain mean basically zero, to apply phi trend filter or not
    AB_edge_data_info.thrhld_phi_trend = thrhld_phi_trend; % QC threshold define what is a phi angle trend, throw out frames not obeying trend
    AB_edge_data_info.thrhld_strain_range = thrhld_strain_range; % QC threshold strain range, to label suspicious
    AB_edge_data_info.thrhld_ratio_stdMeanStr = thrhld_ratio_stdMeanStr; % QC treshold ratio of strain std to mean, to label suspicious
    AB_edge_data_info.A_suspicious = A_suspicious;      % suspicious = scattered over a wide range
    AB_edge_data_info.A_suspicious = phi_suspicious;    % suspicious = large phi range or same number pos/neg angles while not close to zero angle and sufficiently large strain 
    AB_edge_data_info.A_suspicious = strain_suspicious; % suspicious = large std compared to mean and have a wide range.
    
    % Manual output filename entery
%     fprintf('\n%s\n',srcFiles(1).name);    
%     output_filename = input('Output filename:   (AB ID automatically included)    ','s');
    % -- Update 9/10/21: Auto output file name -- %
    % Assumed multipage tif with file name format includes '_MMStack_'
    % Starting location of '_MMStack_' in tiff filename
    mov_name = srcFiles(1).name;
    mm_txt = strfind(mov_name,'_MMStack_');
    % Pickout the text before '_MMStack_'
    output_filename = strcat(qual,'_',mov_name(1:mm_txt-1));
    fprintf('Output f/n: \t%s\n\n',output_filename); 
    % -- Update 9/10/21: Auto output file name -- %
    
    % Save MAT data and shape metric figures. Prefix 'AB<ID>_shape_
    dir_output_filename = strcat(datafolder,'AB',sprintf('%02d',ID),...
        '_shape_',output_filename); % for Mac
    save(dir_output_filename,'AB_edge_data','AB_edge_data_info')
    dir_fig_filename1 = strcat(datafolder,'AB',sprintf('%02d_',ID),output_filename,...
        '_fig_perim_area.png'); % for Mac
    dir_fig_filename2 = strcat(datafolder,'AB',sprintf('%02d_',ID),output_filename,...
        '_fig_ab_phi_strain.png'); % for Mac
    dir_fig_filename5 = strcat(datafolder,'AB',sprintf('%02d_',ID),output_filename,...
        '_fig_edge_ellipse_lastFrame.png'); % for Mac
    saveas(hfig(1), dir_fig_filename1) % area, perimeter scatter plot
    saveas(hfig(2), dir_fig_filename2) % ellipse axes, strain scatter plot
    saveas(hfig5,   dir_fig_filename5) % last plot check of detected contour over stack frame
end


clear A_suspicious phi_suspicious strain_suspicious numFr_suspicious

%% Quick shape measurement
strain_mean = mean(AB_edge_data.strain);
strain_SEM = std(AB_edge_data.strain)/sqrt(length(AB_edge_data.strain));
dia_mean = mean([mean(AB_edge_data.a) mean(AB_edge_data.b)])/9.3457; % um, assumes 100x objective on Zeiss
fprintf('AB %i Quick shape measurement:\n', ID)
fprintf('   Average strain:  %6.4f \n', strain_mean)
fprintf('   SEM strain:      %6.4f \n', strain_SEM)
fprintf('   Average diameter [um]:  %6.4f \n\n\n', dia_mean) 


%% END AB EDGE DETECTION LOOP
end % loop over each AB trajectory

pause(5)
beep

disp('run_apoptoticBody_shape_traj COMPLETE!')

%% Copyright conditions

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.