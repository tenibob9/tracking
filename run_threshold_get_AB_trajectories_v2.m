%% Determine centroids of detected potential L-EVs for stiffness analysis
% 
% V2 8/27/21: For batch processing
% - Automatically determine file name from text preceeding '_MMStack' in
%   tif file name
% - User selects folders to analyze
% - Optional modification of trajectories
%
% Terminology changes Feb 28, 2023
% - Changed display text to 'L-EV' from 'AB'. Variable names were not
% changed.
% AB = apoptotic body, which is too specific. L-EV = large extracellular
% vesicle which is a general term for vesicles that are >1 µm diameter

close all; clear; clc

%% ANALYSIS PARAMETERS
% May need to be optimized for each sample, experiment, and/or
% microscope/camera/objective system

% IMAGE PREPROCESSING PARAMETERS
% standard deviation of Gaussian filter with 2D smoothing kernal
std_gfilt = 2; % px

% THRESHOLDING PARAMETERS
binWidth = 9;   % good bin width for meaningful histogram for 3/26/21 movies
% Assumes 16-bit image
% [counts,binLocations]=imhist(Data_bsub_smth(:,:,100),round((2^16-1)/binWidth));
% For L-EV movies, a global threshold of too lax 280, too many false positives
% - 335 seems to work well for 3/26/21 movies
% - 280 works well for 8/23 and 8/25/21 MGG experiments
% - 280 works for 10/29/21, 11/12/21 Gli36
% - Tried 260, 270 for 2/24/23 BrCa HCC1937 expts, but picked up too much crap. Back to 280
gl_thrhld = 280; 
%gl_thrhld = 300;  % for 1000uLhr Sara

% SHAPE PARAMETERS
% Looking for objects that have an area 0.5 - 10 µm^2 (~ 40 - 875 px^2)
% and circularity >0.9 (for 100x objective on Zeiss Axiovert)
minC = 0.9;
maxC = 1.1;
minA = 40;  % px^2
maxA = 875; % px^2

% PARTICLE TRACKING PARAMETERS
% Particle tracking max displacment: an estimate of the maximum distance 
% that a particle would move in a single time interval.
% •	Maximum displacement (maxdisp) parameter of particles in a trajectory 
% for track.m depends on the flow rate/average velocity of L-EVs
% Smaller maxdisp means less manual fixing of trajectories
% o	Maxdisp = 75 px, even 60 px (best), good for low flow rates < 1250 µL/hr (1 mm wide, 200 µm
%   deep = mean velocity < 2.3 mm/s).
% o	Maxdisp  = 90 px (best) (or 100 px) better for ≥ 2.3 mm/s
% o	Maxdisp = 110 for ≥ 1750 µL/hr (150 px was too big and led to erroneous,
%   unrealistic trajectories)
maxdisp = 75; % px, reduce distance if 'Excessive Combinitorics' error from tracks
%maxdisp = 55; % px, I used 30 for 3/31/23 expt. z=40 movie 1. I tried: 110-60-50-40-35 (there was progress)-then 30 
minFrames = 5;  % minimum number of frames in trajectory
% Distance threshold to combine trajectories
% If particles are closer that this distance in a given frame, combine the
% trajectoies because they are the same L-EV.
% Should be on the order of the average particle size, not much bigger
d_combo = 15;    % px

%% THRESHOLDED experimental image files from movies

%% Preprocessing
% Previously done in ImageJ_AB_thersholding.ijm
% optional plot at intermediary steps
iplot = 0;

%% Load movie frames into 3D array
% Load original experimental images
% Specify directory
disp('Select directory containing the ORIGINAL movie frames:')
directory = strcat(uigetdir('', 'Select directory containing the ORIGINAL movie frames:'),'/');
% Find the movie tiff files in that directory
srcFiles = dir( sprintf('%s//*.tif', directory) );  % look for tiffs
%srcFiles(1) = []; % had to do this for 750uLhr because my modified was
%also saved and I couldn't delete. and 1250uLhr
N = length(srcFiles); % number of tiff files

if N == 1       % multi-page tiff
    flag_multi = 1;     % turn on flag for multipage tiff source files
    mov_info = imfinfo(strcat(directory,srcFiles.name));
    numframes   = length(mov_info);   % number of frames
    mov_data_type = strcat('uint',num2str(mov_info(1).BitDepth));    % 8-bit or 16-bit images
    Data_stack = zeros(mov_info(1).Height, mov_info(1).Width, ...227:228
        numframes, mov_data_type); % preallocate
elseif N > 1    % individual tiff files (image sequence)
    numframes = N;      % number of frames
    flag_multi = 0;     % multipage tiff flag off
    disp('This code only works for multipage tiffs.')
else
    disp('Problem reading in original movie tiff files.')
    return
end

%%  Read images into MATLAB 3D array

for nn = 1:numframes
    if flag_multi
        % Get filename, and load image for multipage tiff movies
        filename = strcat(directory, srcFiles.name); 
        image = imread(filename,nn);    % Load image
        Data_stack(:,:,nn) = image;
    else
        % Get filename, and load image for movies as individual tiff files
        filename = strcat(directory, srcFiles(nn).name);
        image = imread(filename);   % Load image
    end
end

[rows, columns, numSlices] = size(Data_stack);

if iplot, figure, imshow(Data_stack(:,:,round(numSlices/2)),'DisplayRange',[]), title('Checking movie import. This is a frame in the middle of the stack.'), end

%% Z projection - average image = background for this movie
% run("Z Project...", "projection=[Average Intensity]");

bkgrdImg = zeros(rows, columns, class(Data_stack)); % Or whatever class you want.
for col = 1 : columns
    for row = 1 : rows
        thisZVector = Data_stack(row, col, :);
        meanValue = mean(thisZVector);
        bkgrdImg(row, col) = meanValue;
    end
end
disp('Background image found.')

if iplot, figure, imshow(bkgrdImg), title('Background = ZProj Average'), end


%% Subtract background then smooth image with Gaussian filter
% It is very important to do smoothing before thresholding to avoid speckles
% Reproducing these ImageJ macro commands
% imageCalculator("Subtract create stack", "file1.tif","AVG_stack.tif");
% run("Gaussian Blur...", "sigma=2 stack");

% preallocate
Data_bsub      = zeros(size(Data_stack),class(Data_stack));
Data_bsub_smth = zeros(size(Data_stack),class(Data_stack));



% Store axes handles here to fix a plotting timing issue later
if iplot, f = figure; ax = axes(f); imshow(Data_bsub(:,:,round(numSlices/2)),'DisplayRange',[],'Parent',ax), title(ax,'Background-subtracted Frame'), end
if iplot, f2 = figure; ax2 = axes(f2); imshow(Data_bsub_smth(:,:,round(numSlices/2)),'DisplayRange',[],'Parent',ax2), title(ax2,'Background-subtracted, Smoothed Frame'), end

%% Threshold movie slices with a global threshold

% Based on these ImageJ macro commands
% run("Threshold...");  // open thresholding window
% setAutoThreshold("Otsu dark stack");
% call("ij.plugin.frame.ThresholdAdjuster.setMode", "B&W");
% setSlice(112);
% //setThreshold(607, 65535); // code up to this point seemed to make good thresholding choice
% setOption("BlackBackground", false); // doesn't seem to do anything here, for Apply window
% run("Convert to Mask", "method=Otsu background=Dark black list");


% imbinarize function needs a scaled threshold between 0 and 1
if isa(Data_bsub_smth,'uint16')
    T = gl_thrhld/(2^16-1);     % 16-bit images
elseif isa(Data_bsub_smth, 'uint8')
    T = gl_thrhld/(2^8-1);      % 8-bit images
else
    disp('Unsupported image class. Try again')
end

BW = zeros(size(Data_stack),'logical'); % preallocate

% Threshold movie slices to global threshold
for ii = 1:numSlices
    BW(:,:,ii) = imbinarize(Data_bsub_smth(:,:,ii), T);
end
disp('Thresholding complete!')

if iplot, f = figure; ax = axes(f); imshow(BW(:,:,1),'Parent',ax), title(ax,'Thresholded Frame'), end


%% Find centroids of detected particles in a certain size range and circularity

fid = fopen('potential_AB_centroids.txt','w');
pos_list = [];

for nn = 1:numframes

    image = BW(:,:,nn);     % current image

    [B,L] = bwboundaries(image);

    if iplot
        figure(1), imshow(image); hold on
        title(['Frame ',num2str(nn)])
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2),boundary(:,1),'b.')    % imshow(image)
        %     plot(boundary(:,1),boundary(:,2),'b.')      % imshow(image')
        end
    end

    stats = regionprops(L,'Area','Centroid');

    % loop over the boundaries. Estimate object area and perimeter using 2
    % approaches: (1) integrating spline estimates of object shape and 
    % (2) simple
    for k = 1:length(B)

      % obtain (X,Y) boundary coordinates corresponding to label 'k'
      boundary = B{k};

      % skip to next iteration of "boundary" is just two points
      if size(boundary,1) < 3
          continue
      end

        % Use Splines to estimate the object area and perimeter
        % ------ polar coordinates (r(theta),theta) ------ %
        % Convert (x,y) coordinates to (r(theta),theta) coordinates
        x = boundary(:,1);
        y = boundary(:,2);
        [rho, theta, x_new, y_new] = get_rho_theta(x, y);
        if length(rho) < 3
            continue
        end
        % Code from revolved_area_vol_spline.m in vesicle code MATLAB Dropbox
        smth = 0.9;
        scs = csaps(theta, rho, smth);   % smooth curve between scattered data points
        %   scs = csape(theta, rho, 'clamped'); % goes through every point
        spline_pts_polar = fnplt(scs);  % piecewise polynomial [theta rho(theta)]
        theta_spline = spline_pts_polar(1,:);
        rho_spline = spline_pts_polar(2,:);
        x_spline = rho_spline.*cos(theta_spline)+ mean(x_new); % global coord
        y_spline = rho_spline.*sin(theta_spline)+ mean(y_new); % global coord
        if iplot
            figure(1)
            plot(y_spline, x_spline,'-r')
            plot([y_spline y_spline(1)], [x_spline x_spline(1)],'-m')
        end

        % Supply details of the piecewise polynomial:
        %     breaks, coefficients, number of pieces, order and dimension of target
        [breaks,coefs,l,kk,d] = unmkpp(scs);
        % Make the polynomial that describes the derivative drho/dtheta
        % Smooth cubic spline differentiated for drho/dtheta vs dtheta
        scsd = mkpp(breaks,repmat(kk-1:-1:1,d*l,1).*coefs(:,1:kk-1),d);

        % Integrate to determine surface area, volume of revolved curve
        % (2) GAUSSIAN QUADRATURE
        % Set number of Gauss points nq
        nq = 5;
        % Get Gauss points and weights
        [x,w] = gaussian_quad(nq);

        % Initialize
        A = 0; L = 0;

        for ii = 1:l         % loop over piecewise polynomial segments

            % interval end points. Here z = theta, the independent coord.
            z0 = scs.breaks(ii);
            z1 = scs.breaks(ii+1);
            % Jacobian of mapping dz/dx in this interval
            Jq = (z1-z0)/2;

            for q = 1:nq    % loop over Gauss points
                % q = Gauss point index
                % x(q) = Gauss point 
                % w(q) = weight

                % Gauss point in transformed to z domain
                zq = (z1-z0)/2*x(q) + (z1+z0)/2;

                % r(zq(x(q))) radius at Gauss point
                rq = polyval(scs.coefs(ii,:), zq - scs.breaks(ii));
                % drdzq(zq(x(q))) radius at Gauss point
                drdzq = polyval(scsd.coefs(ii,:), zq - scsd.breaks(ii));

                % Integrate by adding on Gauss point
                A = A + 0.5*rq^2*w(q)*Jq;
                L = L + sqrt(rq^2 + drdzq^2)*w(q)*Jq;

            end
        end

      % (2) Simple estimate of the object's perimeter
      % add up length of line segments connecting neighboring points
      delta_sq = diff(boundary).^2;    
      perimeter = sum(sqrt(sum(delta_sq,2)));
      % Simple estimate of area  corresponding to label 'k'
      area = stats(k).Area;

      % Compare area, perimeter measurements from the 2 approaches
      area_diff = (A-area)/area;
      perim_diff = (L-perimeter)/perimeter;

      % compute the roundness metric, circularity
      Cest = 4*pi*area/perimeter^2;     % conclusion: too rough
      C_spline = 4*pi*A/L^2;            % a better metric for L-EVs

      % display the circularity results from both approaches
      metric_string = sprintf('B %i %2.2f %2.2f',k,Cest,C_spline);

      % mark objects above the circularity and area thresholds with a red circle
      % Using spline perimeter, area estimates
      if (C_spline >= minC) && (C_spline <= maxC) && (A >= minA) && (A <= maxA)
        centroid = stats(k).Centroid;
        if iplot
            plot(centroid(1),centroid(2),'ro');
            text(boundary(1,2)-15,boundary(1,1)+5,metric_string,'Color','y',...
               'FontSize',14,'FontWeight','bold')
        end
        % Add this object to the list of potential L-EVs for shape analysis
        fprintf(fid,'%i \t%i \t%4.3f \t%4.3f \t%4.3f \t%4.3f \t%7.4f \t%7.4f\n', ...
            nn, k, area_diff, perim_diff, Cest, C_spline, centroid(1), centroid(2));
        % Save detected particle positions and frame number
        % Define the track function input 'pos_list'. First two columns are x and y
        % locations of each particle. The third column is the frame number of
        % that particle in that location.
        pos_list = [pos_list; centroid(1) centroid(2) nn];
      end

%       % Using rough area, perimeter estimates
%       if (Cest >= minC) && (Cest <= maxC) && (area >= minA) && (area <= maxA)
%         centroid = stats(k).Centroid;
%     %     plot(centroid(1),centroid(2),'c*');
%       end

    end

    % pause
    hold off

end

disp('Centroids of potential L-EVs found. Starting to form trajectories from these centroids...')

fclose(fid);

%% Link particle locations to form trajectories
% http://www.physics.emory.edu/faculty/weeks//idl/

tr = track(pos_list, maxdisp);
% Output is a matrix of 4 columns. The first 3 columns are the x, y and
% frame of the particles sorted into a series of trajectories. The last
% column is the 'id number' which is the identified particle trajectory.
% tr = [x_COM y_COM frame_num particle_num]

%%  Visual output
hfig = figure(2); % all trajectories
plot_traj(tr, image)
xlabel('X [pixel]')
ylabel('Y [pixel]')
title(['Particle Trajectories max disp ',num2str(maxdisp),...
    ' px , ',num2str(tr(end,end)),' trajectories for ',num2str(size(pos_list,1)),...
    ' particles'])

fontSize = 18;
% Set font sizes and line widths on axes and its labels
set(gca,"FontSize",fontSize,"LineWidth",2,"FontWeight","bold","TickLength",[0.025 0.025]);
% Set all other text in figure to fontsize
set(findall(gcf,"type","text"),"FontSize",fontSize);

disp('Trajectories found. Now filtering trajectories...')

%% Filter trajectories
% By number of frames in trajectory

tr_filt = tr;   % new variable so do not overwrite full track list
for ii = 1:tr(end,4)
    particle_tr_numFr(ii) = sum(tr(:,4)==ii);
    if particle_tr_numFr(ii) < minFrames
        tr_filt(tr_filt(:,4)==ii,:) = [];     % get rid of full rows
    end
    particle_trFilt_numFr(ii) = sum(tr_filt(:,4)==ii);
end

hall_stem = figure(3); stem(particle_tr_numFr), hold on
%     plot([tr(1,4) tr(end,4)], 2*[1 1], 'k--')
xlabel('Particle number from tr')
ylabel('Number of frames in trajectory')
title('Raw Trajectories')

hfilt = figure(4); stem(particle_trFilt_numFr), hold on
plot([tr_filt(1,4) tr_filt(end,4)], minFrames*[1 1], 'k--') % minFrame threshold
legend('Trajectories','Min. num frames')
xlabel('Particle number from tr')
ylabel('Number of frames in trajectory')
title('Filtered Trajectories: Min. Frames in a Trajectory')

% Revise trajectory title
figure(hfig)
title(['Particle Trajectories max disp ',num2str(maxdisp),...
' px , ',num2str(length(unique(tr_filt(:,4)))),' filtered trajectories for ',num2str(size(pos_list,1)),...
' particles detected'])

num_filt_traj = length(unique(tr_filt(:,4)));
disp(['Particle Trajectories max disp ',num2str(maxdisp),...
' px , ',num2str(length(unique(tr_filt(:,4)))),' filtered trajectories for ',num2str(size(pos_list,1)),...
' particles detected'])


%% Combine trajectories: automatically based on frame and distance

% tr = [x_COM y_COM frame_num particle_num]

tr_filt_vis = tr_filt;   % new variable so do not overwrite minFrame filtered track list
tr_filt_vis_preAutoCombo = tr_filt_vis; % save copy before combining L-EVs

AB_IDs = unique(tr_filt_vis(:,4));

for ind_ID = 1:length(AB_IDs) % loop over each L-EV trajectory
    
    ID = AB_IDs(ind_ID);    % L-EV particle number
    % Pick out rows in tr_filt that correspond to L-EV with current ID
    AB_traj_info = tr_filt_vis( tr_filt_vis(:,4)==ID , :); 
    x_ID = AB_traj_info(:,1);       % approx. AB x coordinate [px]
    y_ID = AB_traj_info(:,2);       % approx. AB y coordinate [px]
    frames_ID = AB_traj_info(:,3);  % trajectory frames
    
    % preallocate combo ID list with this L-EV
    combo_IDs = ID;
    
    for jj = 1:length(frames_ID) % loop over frames for this trajectory
        
        fr = frames_ID(jj);     % frame in trajectory
        x_fr = x_ID(jj); y_fr = y_ID(jj);   % ID position at frame
        % Identify other particles in this frame
        particles_same_frame = tr_filt_vis( tr_filt_vis(:,3)==fr ,:);      % rows at same frame number
        particles_same_frame( particles_same_frame(:,4)==ID, :) = [];      % throw out rows from current ID trajectory
        % Particle trajectory numbers (should all be unique already).
        % Need stable to keep AB_IDs in same order as in
        % particles_same_frame.
        % Otherwise unintentional combining of wrong L-EVs.
        AB_IDs_same_frame = unique(particles_same_frame(:,4),'stable');    
        x_same_frame = particles_same_frame(:,1);    % x coordinates
        y_same_frame = particles_same_frame(:,2);    % y coordinates
        % check distance to each particle in frame
        for kk = 1:length(x_same_frame)
            % Distance is the 2-norm
            d = norm([x_fr y_fr] - [x_same_frame(kk) y_same_frame(kk)]); % px
            if d < d_combo  % if close enough, then these are the same L-EV
                AB_IDs_same_frame(kk);
                combo_IDs = [combo_IDs; AB_IDs_same_frame(kk)]; % append to combo list
            end
        end
            
    end     % end loop of L-EV ID's frames
    
    combo_IDs = unique(combo_IDs,'stable');
    
    if length(combo_IDs)>1 % copy and paste form combine_del_trajectories
        x_tr = []; y_tr = []; frames_tr = [];
        for ii = 1:length(combo_IDs)
            % Trajectory info for this L-EV
            ID2 = combo_IDs(ii);    % particle number
             % Pick out rows in tr_filt that correspond to L-EV with current ID
            AB_traj_info2 = tr_filt_vis( tr_filt_vis(:,4)==ID2 ,:);
            x_tr = [x_tr; AB_traj_info2(:,1)];       % approx. AB x coordinate [px]
            y_tr = [y_tr; AB_traj_info2(:,2)];       % approx. AB y coordinate [px]
            frames_tr = [frames_tr; AB_traj_info2(:,3)];  % trajectory frames
        end

        % delete duplicate frames, taking first occurance of frame in list
        [C,ia,~]  = unique(frames_tr);   % find unique frames list, outputting indicies ia
        % reduce x_tr, y_tr, frames_tr to unique frames
        x_tr = x_tr(ia);
        y_tr = y_tr(ia);
        frames_tr = frames_tr(ia);

        AB_ID = combo_IDs(1); % Assign ID number of first particle as ID

        % Assemble this combo trajectory in tr_filt_vis format
        combo_AB = [x_tr y_tr frames_tr AB_ID*ones(size(x_tr))];
        
%         % check not loosing rows
%         combo_IDs
%         whos x_tr y_tr frames_tr AB_ID combo_AB
        tr_filt_vis_PreComboII = tr_filt_vis;   % save before combine these L-EVs
        
        % Delete rows of combined L-EVs from tr_filt_vis
        for ii = 1:length(combo_IDs)
            % Trajectory info for this L-EV
            ID3 = combo_IDs(ii);    % particle number
            tr_filt_vis(tr_filt_vis(:,4)==ID3, :) = [];  % get rid of full rows
        end
        
        % Append combined trajectory onto tr_filt_vis
        tr_filt_vis = [tr_filt_vis; combo_AB]; 
        
        % plot check
        figure(20); 
        subplot(1,2,1), plot_traj(tr_filt_vis,image); % BEFORE
        subplot(1,2,2), plot_traj(tr_filt_vis_PreComboII,image); % AFTER

    end     % end combining code

end     % end loop over L-EV trajectories

%% Visual inspection of detected trajectories: Keep or discard

% Code needed if visual inspection is before auto L-EV combining
% tr_filt_vis = tr_filt;   % new variable so do not overwrite minFrame filtered track list
% filt_traj_list = unique(tr_filt(:,4));
% n_filt_traj = length(filt_traj_list);

filt_traj_list = unique(tr_filt_vis(:,4));
n_filt_traj = length(filt_traj_list);

% Update 12/9/21 Save zoom snapshot with filename indicating if kept or discarded
% Save in subfolder AB_tracks_zoom_imgs
[success, message, messageID] = mkdir(directory, 'AB_tracks_zoom_imgs');

for ii = 1:n_filt_traj
    
    % Pick out trajector info for this L-EV
    ID = filt_traj_list(ii);    % particle number
    
    % Pick out rows in tr_filt that correspond to L-EV with current ID
    AB_traj_info = tr_filt_vis( tr_filt_vis(:,4)==ID ,:); 
    x_tr = AB_traj_info(:,1);       % approx. AB x coordinate [px]
    y_tr = AB_traj_info(:,2);       % approx. AB y coordinate [px]
    frames_tr = AB_traj_info(:,3);  % trajectory frames
    
    % Plot original images with trajectory COM points (evenly spaced minFrames
    % frames throughout trajectory - update 12/7/21)
    ind_check = round(linspace(1, length(frames_tr), minFrames));
    % Update 12/9/21: Arrange figures in 2 rows. Tried mod but not clever
    % enough. Just hard code for minFrames = 5
    i_xcorner = [0 1 2 0 1]; i_ycorner = [1 1 1 0 0];
    for jj = 1:minFrames
        % Load current frame jj
        if flag_multi
            % Get filename, and load image for multipage tiff movies
            filename = strcat(directory, srcFiles.name); 
%             image = imread(filename,frames_tr(jj));    % Load image, old first minFrames frames
            image = imread(filename,frames_tr(ind_check(jj)));    % Load image, 12/6/21 update evenly spaced
        else
            % Get filename, and load image for movies as individual tiff files
            filename = strcat(directory, srcFiles(frames_tr(ind_check(jj))).name); % 12/6/21 update evenly spaced
%             filename = strcat(directory, srcFiles(frames_tr(jj)).name); % old first minFrames frames
            image = imread(filename);   % Load image
        end
        hfig_check = figure(10+jj); ax_check = gca; % Assign higher figure numbers
%         set(gcf,'Position',[(1+500*(jj-1)) 680 525 525]) % one row of figures
        %set(gcf,'Position',[(1+500*i_xcorner(jj)) (1+680*i_ycorner(jj)) 525 525]) % update 12/9/21 2 rows of figures
        %set(gcf,'Position',[(1+500*i_xcorner(jj)) (1+680*i_ycorner(jj)) 380 380]) % Sara's mac window
        set(gcf,'Position',[(1+500*i_xcorner(jj)) (1+500*i_ycorner(jj)) 525 525]) % Sara's mac window bigger view
        imshow(image,[],'Parent',ax_check);
%         % old: consecutive frames 1 - minFrames
%         hold on; plot(x_tr(jj),y_tr(jj),'yo','MarkerSize',20);
%         title(ax_check, ['AB Trajectory #', num2str(filt_traj_list(ii)), ...
%             ', Global Frame ',num2str(frames_tr(jj))]); 
        % 12/7/21 update: evenly spaced across trajectory
        hold on; plot(x_tr(ind_check(jj)),y_tr(ind_check(jj)),'yo','MarkerSize',20);
        title(ax_check, ['L-EV Trajectory #', num2str(filt_traj_list(ii)), ...
            ', Global Frame ',num2str(frames_tr(ind_check(jj)))]); 
        hold off
        
        % at midpoint, plot a zoomed in image
        if jj == round(minFrames/2)  % at midpoint, plot a zoomed in image
            hfig_zoom = figure(20+jj); ax_check2 = gca; % Assign higher figure numbers
            %set(gcf,'Position',[(1+500*(jj-1)) 15 525 525])
            %set(gcf,'Position',[(1+500*(jj-1)) 15 380 380])   % Sara's mac window
            set(gcf,'Position',[(1+500*(jj-1)) 15 525 525])   % Sara's mac window bigger image
            imshow(image,[],'Parent',ax_check2);
%             imshow(image,[]);
            % Zoom in by setting axis limits to 10% of image size, centered
            % at the particle center point
            Lx = size(image,1)/10; Ly = size(image,2)/10;
%             % old: consecutive frames 1 - minFrames
%             xlim([x_tr(jj)-Lx x_tr(jj)+Lx]);
%             ylim([y_tr(jj)-Ly y_tr(jj)+Ly]);
%             title(ax_check2, ['ZOOM: AB Trajectory #', num2str(filt_traj_list(ii)), ...
%             ', Global Frame ',num2str(frames_tr(jj))]);
            % 12/7/21 update: evenly spaced across trajectory
            xlim([x_tr(ind_check(jj))-Lx x_tr(ind_check(jj))+Lx]);
            ylim([y_tr(ind_check(jj))-Ly y_tr(ind_check(jj))+Ly]);
            title(ax_check2, ['ZOOM: L-EV Trajectory #', num2str(filt_traj_list(ii)), ...
            ', Global Frame ',num2str(frames_tr(ind_check(jj)))]);
            
        end
    end
    
    flag_keep = input('Keep this L-EV trajectory? 1 = yes, 0 = no      ');
        
    if flag_keep
        fprintf('Keeping L-EV #%i (%i out of %i)\n',filt_traj_list(ii),ii,n_filt_traj); % update 12/9/21 display more info
        % Update 12/9/21 Save zoom snapshot with filename indicating if kept or discarded
        dir_fig_filename_zoom = strcat(directory,'AB_tracks_zoom_imgs/','fig_zoom_KEEP_AB',num2str(filt_traj_list(ii)),'.png'); % for Mac
    else
        fprintf('Discarding L-EV #%i (%i out of %i)\n', filt_traj_list(ii),ii,n_filt_traj); % update 12/9/21 display more info
        tr_filt_vis(tr_filt_vis(:,4)==ID, :) = [];  % get rid of full rows
        % Update 12/9/21 Save zoom snapshot with filename indicating if kept or discarded
        dir_fig_filename_zoom = strcat(directory,'AB_tracks_zoom_imgs/','fig_zoom_DISCARD_AB',num2str(filt_traj_list(ii)),'.png'); % for Mac
    end
    % Update 12/9/21 Save zoom snapshot with filename indicating if kept or discarded
    saveas(hfig_zoom, dir_fig_filename_zoom)
    
end


%% Plot of automatically filtered trajectories

hfig_filt = figure(5);  pause(0.1)   % trajectories of filtered by frame number
plot_traj(tr_filt,    image)
title('Filtered by number of frames')

% Revise trajectory plot and title
hfig_vis = figure(6); pause(0.1) % trajectories of filtered + visually inspected
plot_traj(tr_filt_vis,image)
title(['Particle Trajectories max disp ',num2str(maxdisp),...
' px , ',num2str(length(unique(tr_filt_vis(:,4)))),' filtered trajectories (minFrames, visual inspection) for ',num2str(size(pos_list,1)),...
' particles detected'])

disp(['Particle Trajectories max disp ',num2str(maxdisp),...
' px , ',num2str(length(unique(tr_filt_vis(:,4)))),' filtered trajectories (minFrames, visual inspection) for ',num2str(size(pos_list,1)),...
' particles detected'])

%% Save detected trajectories

fprintf('\nMovie f/n: \t%s\n',srcFiles(1).name);     
% output_filename = input('Enter a descriptive output filename (ABtracks_ will preceed): ','s');
% --- Update 8/27/21 -- Autoread flow rate from filename --- %
% Assumed multipage tif with file name format includes '_MMStack_'
% Starting location of '_MMStack_' in tiff filename
mov_name = srcFiles(1).name;
mm_txt = strfind(mov_name,'_MMStack_');
% Pickout the text before '_MMStack_'
output_filename = mov_name(1:mm_txt-1);
fprintf('Output f/n: \t%s\n\n',output_filename); 
% --- Update 8/27/21 -- Autoread flow rate from filename --- %

dir_output_filename = strcat(directory, 'ABtracks_', output_filename); % for Mac
save(dir_output_filename,'tr','tr_filt','tr_filt_vis','pos_list','minC',...
    'maxC','minA','maxA','maxdisp','T','gl_thrhld','mov_info','bkgrdImg','d_combo','minFrames')    
[success, message, messageID] = copyfile('potential_AB_centroids.txt', directory);
dir_fig_filename1 = strcat(directory, 'fig_filt_stem.png'); % for Mac
dir_fig_filename2 = strcat(directory, 'fig_filt_stem.fig'); % for Mac
saveas(hfilt, dir_fig_filename1)
saveas(hfilt, dir_fig_filename2)
dir_fig_filename3 = strcat(directory, 'fig_filt_traj.png'); % for Mac
dir_fig_filename4 = strcat(directory, 'fig_filt_traj.fig'); % for Mac
saveas(hfig_vis, dir_fig_filename3)
saveas(hfig_vis, dir_fig_filename4)

disp('L-EV particle tracking complete!')

%% OPTIONAL Modify L-EV trajectories

% STEP 1. Save copy of tr_filt_vis
tr_filt_vis_preComboDel = tr_filt_vis;

% choice_mod = input('Do you need to modify any L-EV trajectories? y = yes or n = no  ','s');
choice_mod = 'y'; % hard code 'yes'

if strcmp(choice_mod,'n') % no modifications needed 
    
    % STEP 3. Update tracks mat file and save bck up
    % Just save backup copy of trajectory data
    % Backup copy
    save(strcat(directory,'BACKUP_tr_filt_vis.mat'),...
        'dir_output_filename','tr_filt_vis','tr_filt_vis_preComboDel') 
    % Append to ABtracks mat file
    save(dir_output_filename,'tr_filt_vis','tr_filt_vis_preComboDel','-append') 
    disp('Modified L-EV trajectories saved!')

elseif strcmp(choice_mod, 'y') % Modifications necessary
    
    % STEP 2. Custom Cursor
    % Activate the custom cursor and explore the trajectories to see which ones
    % should be edited using the appropriate section.
    figure(hfig_vis)
    dcm_obj = datacursormode(hfig_vis);
    dcm_obj.Enable = 'on';
    dcm_obj.UpdateFcn = {@myupdatefcn,tr_filt_vis};
    % set(dcm_obj,'UpdateFcn',{@myupdatefcn,tr_filt_vis})
    
    % Modification options
    fprintf('\nL-EV track modification options:\n');
    fprintf('\tf = Delete FRAMES from a trajectory\n');
    fprintf('\tc = COMBINE L-EV trajectories\n');
    fprintf('\tt = Delete an L-EV particle TRAJECTORY\n');
    fprintf('\tq = Quit trajectory modification menu. Done.\n\n');
    task_opt = input('Pick an L-EV track modification option:   ','s');
    
    % While loop
    while ~strcmp(task_opt,'q')     % while not choosing to exit
        
        if strcmp(task_opt,'f')     % Delete FRAMES from a trajectory
            
            %% AS NEEDED: Delete FRAMES from a trajectory
            % Run this section for every trajectory that needs frames to be deleted
            delFr_ID = input('Enter particle ID of trajectory that needs frames removed:   ');

            delFr = input('Enter frame number(s) that needed to be deleted as a vector:   ');

            for ii = 1:length(delFr)
                fr = delFr(ii); % frame to delete
                logDelFr = (tr_filt_vis(:,3)==fr) + (tr_filt_vis(:,4)==delFr_ID);
                [~, im] = max(logDelFr); % find row for this particle ID and this frame
                tr_filt_vis(im, :) = [];  % get rid of full rows
            end

            % Revise trajectory plot and title
            hfig_vis = figure;
            plot_traj(tr_filt_vis,image)
            title(['Particle Trajectories max disp ',num2str(maxdisp),...
            ' px , ',num2str(length(unique(tr_filt_vis(:,4)))),' filtered trajectories (minFrames, visual inspection, combined) for ',num2str(size(pos_list,1)),...
            ' particles detected'])
            dir_fig_filename1 = strcat(directory, 'fig_trajectories_filt_vis_combo.png'); % for Mac
            dir_fig_filename2 = strcat(directory, 'fig_trajectories_filt_vis_combo.fig'); % for Mac
            saveas(hfig_vis, dir_fig_filename1)
            saveas(hfig_vis, dir_fig_filename2)
            
            % Custom cursor
            figure(hfig_vis)
            dcm_obj = datacursormode(hfig_vis);
            dcm_obj.Enable = 'on';
            dcm_obj.UpdateFcn = {@myupdatefcn,tr_filt_vis};
            
            % Modification options
            fprintf('\nL-EV track modification options:\n');
            fprintf('\tf = Delete FRAMES from a trajectory\n');
            fprintf('\tc = COMBINE L-EV trajectories\n');
            fprintf('\tt = Delete an L-EV particle TRAJECTORY\n');
            fprintf('\tq = Quit trajectory modification menu. Done.\n\n');
            task_opt = input('Pick an L-EV track modification option:   ','s');
            
        elseif strcmp(task_opt,'c')     % COMBINE L-EV trajectories
            
            %% AS NEEDED: COMBINE L-EV trajectory sets
            % Run this section for each set of particle trajectories to combine

            combo_IDs = input('Enter the particle IDs to combine as a vector:  ');

            x_tr = []; y_tr = []; frames_tr = [];

            for ii = 1:length(combo_IDs)
                % Trajectory info for this L-EV
                ID = combo_IDs(ii);    % particle number
                 % Pick out rows in tr_filt that correspond to L-EV with current ID
                AB_traj_info = tr_filt_vis( tr_filt_vis(:,4)==ID ,:);
                x_tr = [x_tr; AB_traj_info(:,1)];       % approx. AB x coordinate [px]
                y_tr = [y_tr; AB_traj_info(:,2)];       % approx. AB y coordinate [px]
                frames_tr = [frames_tr; AB_traj_info(:,3)];  % trajectory frames
            end

            % delete duplicate frames, taking first occurance of frame in list
            [C,ia,~]  = unique(frames_tr);   % find unique frames list, outputting indicies ia
            % reduce x_tr, y_tr, frames_tr to unique frames
            x_tr = x_tr(ia);
            y_tr = y_tr(ia);
            frames_tr = frames_tr(ia);

            AB_ID = combo_IDs(1); % Assign ID number of first particle as ID

            % Assemble this combo trajectory in tr_filt_vis format
            combo_AB = [x_tr y_tr frames_tr AB_ID*ones(size(x_tr))];

            % Delete rows of combined ABs from tr_filt_vis
            for ii = 1:length(combo_IDs)
                % Trajectory info for this AB
                ID = combo_IDs(ii);    % particle number
                tr_filt_vis(tr_filt_vis(:,4)==ID, :) = [];  % get rid of full rows
            end

            % Append combined trajectory onto tr_filt_vis
            tr_filt_vis = [tr_filt_vis; combo_AB];

            % Revise trajectory plot and title
            hfig_vis = figure;
            plot_traj(tr_filt_vis,image)
            title(['Particle Trajectories max disp ',num2str(maxdisp),...
            ' px , ',num2str(length(unique(tr_filt_vis(:,4)))),' filtered trajectories (minFrames, visual inspection, combined) for ',num2str(size(pos_list,1)),...
            ' particles detected'])
            dir_fig_filename1 = strcat(directory, 'fig_trajectories_filt_vis_combo.png'); % for Mac
            dir_fig_filename2 = strcat(directory, 'fig_trajectories_filt_vis_combo.fig'); % for Mac
            saveas(hfig_vis, dir_fig_filename1)
            saveas(hfig_vis, dir_fig_filename2)
            
            % Custom Cursor
            figure(hfig_vis)
            dcm_obj = datacursormode(hfig_vis);
            dcm_obj.Enable = 'on';
            dcm_obj.UpdateFcn = {@myupdatefcn,tr_filt_vis};
            
            % Modification options
            fprintf('\nL-EV track modification options:\n');
            fprintf('\tf = Delete FRAMES from a trajectory\n');
            fprintf('\tc = COMBINE L-EV trajectories\n');
            fprintf('\tt = Delete an L-EV particle TRAJECTORY\n');
            fprintf('\tq = Quit trajectory modification menu. Done.\n\n');
            task_opt = input('Pick an L-EV track modification option:   ','s');
            
        elseif strcmp(task_opt,'t')     % Delete an L-EV particle TRAJECTORY
            
            %% AS NEEDED: Delete L-EV PARTICLES because a duplicate or crap
            % Only have to run once if have list of ID of particles to delete
            del_IDs = input('Enter particle IDs to delete as a vector:   ');

            for ii = 1:length(del_IDs)
                ID = del_IDs(ii); % particle ID number
                tr_filt_vis(tr_filt_vis(:,4)==ID, :) = [];  % get rid of full rows
            end

            % Revise trajectory plot and title
            hfig_vis = figure;
            plot_traj(tr_filt_vis,image)
            title(['Particle Trajectories max disp ',num2str(maxdisp),...
            ' px , ',num2str(length(unique(tr_filt_vis(:,4)))),' filtered trajectories (minFrames, visual inspection, combined) for ',num2str(size(pos_list,1)),...
            ' particles detected'])
            dir_fig_filename1 = strcat(directory, 'fig_trajectories_filt_vis_combo.png'); % for Mac
            dir_fig_filename2 = strcat(directory, 'fig_trajectories_filt_vis_combo.fig'); % for Mac
            saveas(hfig_vis, dir_fig_filename1)
            saveas(hfig_vis, dir_fig_filename2)

            % Custom Cursor
            figure(hfig_vis)
            dcm_obj = datacursormode(hfig_vis);
            dcm_obj.Enable = 'on';
            dcm_obj.UpdateFcn = {@myupdatefcn,tr_filt_vis};
            
            % Modification options
            fprintf('\nL-EV track modification options:\n');
            fprintf('\tf = Delete FRAMES from a trajectory\n');
            fprintf('\tc = COMBINE L-EV trajectories\n');
            fprintf('\tt = Delete an L-EV particle TRAJECTORY\n');
            fprintf('\tq = Quit trajectory modification menu. Done.\n\n');
            task_opt = input('Pick an L-EV track modification option:   ','s');

        else
            
            task_opt = input('Invalid L-EV track modification option. Pick again:   ','s');
        
        end
        
    end    
    
    % STEP 3. Update tracks mat file and save backup
    % Backup copy
    save(strcat(directory,'BACKUP_tr_filt_vis.mat'),...
        'dir_output_filename','tr_filt_vis','tr_filt_vis_preComboDel') 
    % Append to ABtracks mat file
    save(dir_output_filename,'tr_filt_vis','tr_filt_vis_preComboDel','-append') 
    disp('Modified L-EV trajectories saved!')
    
else
    
    disp('Modification choice invalid. No modifications for trajectories in this movie.')

end