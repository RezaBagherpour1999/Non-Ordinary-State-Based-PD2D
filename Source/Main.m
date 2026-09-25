% ========================================================================
% Copyright (c) 2022 by Oak Ridge National Laboratory
% Modifications Copyright (c) 2026 by Reza Bagherpour
% All rights reserved.
%
% This file is part of NOSB-PD2D (extended from PDMATLAB2D).
%
% Distributed under a BSD 3-clause license. For the licensing terms see 
% the LICENSE file in the top-level directory.
%
% SPDX-License-Identifier: BSD-3-Clause
% ========================================================================

% ========================================================================
% Main script for running a PDMATLAB2D simulation
% ========================================================================

% Check if simulation output directory exists and create it otherwise
if ~exist(['../Outputs/' InputDeck], 'dir')
    mkdir(['../Outputs/' InputDeck])
end

% ------------------------------------------------------------------------
%                       Create video file(s)
% ------------------------------------------------------------------------

if flag_DynamicPlotting == 1

    % Check if video flag is defined
    if exist('flag_video','var')

        if flag_video == 1
            % Check if video output directory exists and create it otherwise
            if ~exist(['../Videos/' InputDeck], 'dir')
                mkdir(['../Videos/' InputDeck])
            end

            % Find number of plots: a video is generated for each plot
            [s1,~] = size(PlotSettings);

            % Initialize VideoWriter array
            vidfile = VideoWriter.empty(s1, 0);

            % Loop over plots
            for nplot = 1:s1
                % Read plot field name
                field_name = PlotSettings{nplot,1};

                % Video file name
                video_filename = ['../Videos/' InputDeck '/' field_name '.mp4'];

                % Open video file
                vidfile(nplot) = VideoWriter(video_filename,'MPEG-4');
                vidfile(nplot).FrameRate = video_frate;
                open(vidfile(nplot));
            end
        end

    else

        % Define null video flag
        flag_video = 0;

    end

end

% ------------------------------------------------------------------------
%                    Generate grid and neighbor list
% ------------------------------------------------------------------------

% Check if GridFile variable to load grid is defined
if exist('GridFile','var')

    % Check if grid file exists
    if exist(GridFile,'file') == 2

        % ---------------------------------
        %    Load grid and neighbor list
        % ---------------------------------
        tic
        load(GridFile);
        fprintf('Load grid and neighbor list ..................... = %f (sec) \n',toc)

        % Flag for Rectangular Domain Uniform Grid (RDUG)
        flag_RDUG = 0;

        % ---------------------------------
        %        Check grid inputs
        % ---------------------------------
        % Check if variables exist
        gridvars = ["xx","yy","u_NA","IF_NA","V_NA","r_hat_NA","x_hat_NA","y_hat_NA"];

        for n = 1:length(gridvars)           
            if exist(gridvars(n),'var') == 0
                cd ..
                error('Loaded grid variable %s does not exist.', gridvars(n))
            end
        end

        % Check consistency of xx and yy array dimensions
        if isequal(size(xx),size(yy))

        else
            cd ..
            error('xx and yy arrays have inconsistent dimensions.')
        end

        % Check consistency of u_NA, IF_NA, V_NA, r_hat_NA, x_hat_NA, and y_hat_NA array dimensions
        if isequal(size(u_NA),size(IF_NA),size(V_NA),size(r_hat_NA),size(x_hat_NA),size(y_hat_NA)) 

        else  
            cd ..
            error('u_NA, IF_NA, V_NA, r_hat_NA, x_hat_NA, and y_hat_NA arrays do not have all consistent dimensions.')        
        end

    else

        cd ..
        error('Grid file does not exist.')    
    
    end

else
    
    % ---------------------------------
    %        Generate grid
    % ---------------------------------

    tic
    [xx,yy,x,y,dx,dy,VV,xx1,yy1,M] = GridGenerator(Xo,Xn,Yo,Yn,Nx,Ny,PG);
    fprintf('Generate grid ................................... = %f (sec) \n',toc)

    % Tolerance
    tol = 1E-15;

    % Flag for Rectangular Domain Uniform Grid (RDUG)
    if abs(dx - dy) < tol 
        flag_RDUG = 1;
    else
        flag_RDUG = 0;
    end

    % ---------------------------------
    %     Generate neighbor list
    % ---------------------------------

    tic
    [u_NA,IF_NA,V_NA,r_hat_NA,x_hat_NA,y_hat_NA] = NeighborList(Nx,Ny,xx,yy,xx1,yy1,M,del,dx,dy,VV,omega,AlgName,flag_RDUG);
    fprintf('Generate neighbor list .......................... = %f (sec)\n',toc)

end

% ------------------------------------------------------------------------
%                        Create no-fail mask
% ------------------------------------------------------------------------

if exist('nofailfunc','var')
    mask_nofail = nofailfunc(xx,yy);
else
    mask_nofail = 0*xx;
end

% ------------------------------------------------------------------------
%                       Create pre-notch(es)
% ------------------------------------------------------------------------

% Check if array of pre-notches coordinates is defined
if exist('PreNotchCoordinates','var')

    tic

    % Find number of pre-notches
    [s1,~] = size(PreNotchCoordinates);

    % Loop over pre-notches
    for n = 1:s1
        % Coordinates of one endpoint of the pre-notch
        Xc1 = PreNotchCoordinates(n,1);
        Yc1 = PreNotchCoordinates(n,2);

        % Coordinates of the other endpoint of the pre-notch
        Xc2 = PreNotchCoordinates(n,3);
        Yc2 = PreNotchCoordinates(n,4);

        % Create pre-notch
        [u_NA] = PreNotch(xx,yy,u_NA,Xc1,Yc1,Xc2,Yc2);
    end

    fprintf('Create pre-notch(es) ............................ = %f (sec)\n',toc)

end

% ------------------------------------------------------------------------
%                   Compute peridynamic constants
% ------------------------------------------------------------------------

tic

[c,so] = PDBondConstants(omega,del,E,Go,model,PlanarModel);

fprintf('Compute PD constants ............................ = %f (sec)\n',toc)

% ------------------------------------------------------------------------
%                      Impose initial conditions
% ------------------------------------------------------------------------

tic

% Compute initial displacement for all nodes
v = vofunc(xx,yy);   % x-component of initial displacement
w = wofunc(xx,yy);   % y-component of initial displacement

% Compute initial velocity for all nodes
Vv = Vvofunc(xx,yy); % x-component of initial velocity
Vw = Vwofunc(xx,yy); % y-component of initial velocity

fprintf('Impose initial conditions ....................... = %f (sec)\n',toc)

% ------------------------------------------------------------------------
% Compute initial internal force density and macroelastic energy density
% ------------------------------------------------------------------------

tic

[Fv, Fw, W, sigma_bond_stretch_ui] = ForceEnergyDensity(xx, yy, v, w, E, nu, u_NA, del, IF_NA, V_NA, x_hat_NA, y_hat_NA, flag_RDUG);
fprintf('Compute initial internal force & energy densities = %f (sec)\n',toc)

% ------------------------------------------------------------------------
%                Compute initial body force density
% ------------------------------------------------------------------------

tic

bv = bvfunc(xx,yy,0); % x-component of initial body force density
bw = bwfunc(xx,yy,0); % y-component of initial body force density

fprintf('Compute initial body force density .............. = %f (sec)\n',toc)

% ------------------------------------------------------------------------
%                 Compute denominator of damage ratio
% ------------------------------------------------------------------------

% Check if bond breaking is enabled
if flag_BB == 0

elseif flag_BB == 1
    
    if exist('flag_DamagedPrenotches','var')

        % Compute damage ratio denominator before application of pre-notches
        if flag_DamagedPrenotches == 1

            phiD = sum(V_NA,2);
        
        % Compute damage ratio denominator after application of pre-notches
        elseif flag_DamagedPrenotches == 0
        
            phiD = sum((u_NA>0).*V_NA,2);
        
        else
        
            error('flag_DamagedPrenotches should be 0 or 1.')
        
        end

    else

        % Compute default damage ratio denominator (after application of pre-notches)
        phiD = sum((u_NA>0).*V_NA,2);
    
    end

else

    error('flag_BB should be 0 or 1.')

end

% ------------------------------------------------------------------------
%                         Time integration loop
% ------------------------------------------------------------------------

% Create time vector
tVec = Ti:dt:Tf;

% Number of time steps
Nt = length(tVec);

% Initial time
t = Ti;

tic

% ==========================================
% Node Indices Definitions (Physically Correct)
% In this mesh, lower indices correspond to +y (Top)
% ==========================================
col_Mid         = round(Nx / 2);        % Mid-plane column along x
row_Above       = round(Ny / 4);        % Geometric center row of the upper half-rectangle
row_Center      = round(Ny / 2);        % Central row
row_Below       = round(3 * Ny / 4);    % Geometric center row of the lower half-rectangle

idx_TopEdge     = col_Mid;                               % Top edge center (+y)
idx_AboveCenter = (row_Above  - 1) * Nx + col_Mid;       % Center of upper half-rectangle
idx_Center      = (row_Center - 1) * Nx + col_Mid;       % Center of domain
idx_BelowCenter = (row_Below  - 1) * Nx + col_Mid;       % Center of lower half-rectangle
idx_BottomEdge  = (Ny - 1) * Nx + col_Mid;               % Bottom edge center (-y)
idx_Node1       = 1;                                     % Corner node (Node 1)

% ==========================================
% Preallocation for Speed
% ==========================================
num_steps   = Nt - 1;
num_saved20 = floor(num_steps / 20); % Number of sampling steps (every 20 steps)

% Arrays saved every 20 steps
time_history = zeros(1, num_saved20);

% Y-Displacements (w)
y_TopEdge_history     = zeros(1, num_saved20);
y_AboveCenter_history = zeros(1, num_saved20);
y_Center_history      = zeros(1, num_saved20);
y_BelowCenter_history = zeros(1, num_saved20);
y_BottomEdge_history  = zeros(1, num_saved20);

% Kinetic Energy Density
kinetic_energy_TopEdge    = zeros(1, num_saved20);
kinetic_energy_Center     = zeros(1, num_saved20);
kinetic_energy_BottomEdge = zeros(1, num_saved20);

% Y-Force Densities (Fw)
force_TopEdge_y_history     = zeros(1, num_saved20);
force_AboveCenter_y_history = zeros(1, num_saved20);
force_Center_y_history      = zeros(1, num_saved20);
force_BelowCenter_y_history = zeros(1, num_saved20);
force_BottomEdge_y_history  = zeros(1, num_saved20); 
        
% Strain Energy Density (W)
Strain_energy_TopEdge     = zeros(1, num_saved20);
Strain_energy_AboveCenter = zeros(1, num_saved20);
Strain_energy_Center      = zeros(1, num_saved20);
Strain_energy_BelowCenter = zeros(1, num_saved20);
Strain_energy_BottomEdge  = zeros(1, num_saved20);

% Arrays saved every single time step (Node 1)
w_history_Node1     = zeros(1, Nt-1);
force_history_Node1 = zeros(1, Nt-1);
Time                = zeros(1, Nt-1);

% ==========================================
% Time Integration Loop
% ==========================================
for n = 1:Nt-1

    [v,w,Vv,Vw,Fv,Fw,bv,bw,W,u_NA,sigma_bond_stretch_ui] = TimeIntegrator(c_damp, sigma_bond_stretch_ui, TimeScheme, xx, yy, v, w, Vv, Vw, E, nu, del, Fv, Fw, bv, bw, t, bvfunc, bwfunc, dt, u_NA, IF_NA, V_NA, x_hat_NA, y_hat_NA, rho, flag_RDUG, mask_nofail, flag_BB);
    
    % Display progress in command window
    if mod(n,TimeStepDisplayFrequency) == 0
        fprintf('Time integration (n = %4g / %g ) = %f (sec)\n', n, Nt-1, toc);
    end
    
    % Track Node 1 at every step
    w_history_Node1(n)     = w(idx_Node1);
    force_history_Node1(n) = Fw(idx_Node1);
    Time(n)                = t;

    % Break if numerical instability (NaN) occurs
    if isnan(w(1))
        break;
    end
    
    % Sample history data every 20 steps
    if mod(n, 20) == 0
        
        idx = n / 20; % Direct index mapping
        
        % Strain Energy Density
        Strain_energy_TopEdge(idx)     = W(idx_TopEdge);
        Strain_energy_AboveCenter(idx) = W(idx_AboveCenter);
        Strain_energy_Center(idx)      = W(idx_Center);
        Strain_energy_BelowCenter(idx) = W(idx_BelowCenter);
        Strain_energy_BottomEdge(idx)  = W(idx_BottomEdge);
        
        % Force Density in Y direction
        force_TopEdge_y_history(idx)     = Fw(idx_TopEdge);
        force_AboveCenter_y_history(idx) = Fw(idx_AboveCenter);
        force_Center_y_history(idx)      = Fw(idx_Center);
        force_BelowCenter_y_history(idx) = Fw(idx_BelowCenter);
        force_BottomEdge_y_history(idx)  = Fw(idx_BottomEdge); 
        
        % Displacement in Y direction (w)
        y_TopEdge_history(idx)     = w(idx_TopEdge);
        y_AboveCenter_history(idx) = w(idx_AboveCenter);
        y_Center_history(idx)      = w(idx_Center);
        y_BelowCenter_history(idx) = w(idx_BelowCenter);
        y_BottomEdge_history(idx)  = w(idx_BottomEdge);   
             
        % Kinetic Energy Density (0.5 * rho * Vw^2)
        kinetic_energy_TopEdge(idx)    = 0.5 * rho * Vw(idx_TopEdge).^2;
        kinetic_energy_Center(idx)     = 0.5 * rho * Vw(idx_Center).^2;
        kinetic_energy_BottomEdge(idx) = 0.5 * rho * Vw(idx_BottomEdge).^2;    
        
        % Current simulation time
        time_history(idx) = t;
        
    end

    % Update time
    t = tVec(n+1);

    % ==========================================
    % Field Plotting (Dynamic Plotting)
    % ==========================================
    if flag_DynamicPlotting == 0
        % No dynamic plotting
    elseif flag_DynamicPlotting == 1
        if n == 1 || mod(n-1,DynamicPlotFrequency) == 0
            [s1,~] = size(PlotSettings);
            for nplot = 1:s1
                nfig = nplot;                       
                cnodes_name = PlotSettings{nplot,1};    

                if flag_BB == 0
                elseif flag_BB == 1
                    if strcmp(cnodes_name,'Damage')
                        phi = 1 - sum((u_NA>0).*V_NA,2)./phiD;
                    end
                else
                    error('flag_BB should be 0 or 1.')
                end

                cnodes        = eval(PlotSettings{nplot,2}); 
                ctitle        = PlotSettings{nplot,3};       
                psize         = PlotSettings{nplot,4};       
                climits       = PlotSettings{nplot,5};       
                cmap          = PlotSettings{nplot,6};       
                box           = PlotSettings{nplot,7};       
                configuration = PlotSettings{nplot,8};       

                if strcmp(configuration,'Reference')
                    PlotField(nfig,xx,yy,cnodes,ctitle,psize,climits,cmap,box);
                elseif strcmp(configuration,'Current')
                    PlotField(nfig,xx+v,yy+w,cnodes,ctitle,psize,climits,cmap,box);
                else
                    error('Invalid configuration.');
                end
            end
        end
    else
        error('flag_DynamicPlotting should be 0 or 1.');
    end
end

% ==========================================
% Post-Processing / Final Plots
% ==========================================

% Display verification table for array lengths
lens = [ length(time_history);
         length(force_BottomEdge_y_history);
         length(y_BottomEdge_history);
         length(kinetic_energy_BottomEdge) ];

disp(table(["time", "Fy", "y", "KE"]', lens, ...
    'VariableNames', {'Variable','Length'}));

% 2. Kinetic Energy Density vs Time
figure;
plot(time_history, kinetic_energy_TopEdge, 'b.-', ...
     time_history, kinetic_energy_Center, 'r.-', ...
     time_history, kinetic_energy_BottomEdge, 'k.-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Kinetic Energy Density [J/m^3]');
title('Kinetic Energy Density vs Time');
legend('Top Edge', 'Center', 'Bottom Edge', 'Location', 'best');
grid on;


% 4. Force Density (Y) vs Time
figure;
plot(time_history, force_TopEdge_y_history, 'b.-', ...
     time_history, force_AboveCenter_y_history, 'g.-', ...
     time_history, force_Center_y_history, 'r.-', ...
     time_history, force_BelowCenter_y_history, 'm.-', ...
     time_history, force_BottomEdge_y_history, 'k.-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Force Density y [N/m^3]');
title('Force Density (y) on Nodes over Time');
legend('Top Edge', 'Above Center', 'Center', 'Below Center', 'Bottom Edge', 'Location', 'best');
grid on;


% 6. Displacement (Y) vs Time
figure;
plot(time_history, y_TopEdge_history, 'b.-', ...
     time_history, y_AboveCenter_history, 'g.-', ...
     time_history, y_Center_history, 'r.-', ...
     time_history, y_BelowCenter_history, 'm.-', ...
     time_history, y_BottomEdge_history, 'k.-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Displacement in y [m]');
title('Displacement in y on Nodes over Time');
legend('Top Edge', 'Above Center', 'Center', 'Below Center', 'Bottom Edge', 'Location', 'best');
grid on;

% 7. Node 1 Displacement (w) vs Time
figure;
plot(Time, w_history_Node1, 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Displacement w [m]');
title('Displacement w of Node 1 Over Time');
grid on;

% 8. Node 1 Force Density (w) vs Time
figure;
plot(Time, force_history_Node1, 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Force Density w [N/m^3]');
title('Force Density w of Node 1 Over Time');
grid on;
drawnow;

% ==========================================
% 9. Node Position Topology Plot (Corrected Placement)
% ==========================================
figure;
% Plot all nodes
scatter(xx, yy, 15, 'b', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.4); 
hold on;

tracked_indices = [idx_TopEdge, idx_AboveCenter, idx_Center, ...
                   idx_BelowCenter, idx_BottomEdge, idx_Node1];
               
% Highlight tracked nodes
scatter(xx(tracked_indices), yy(tracked_indices), 45, 'r', 'filled', 'MarkerEdgeColor', 'k');

% Verified coordinate-based text placement
text(xx(idx_TopEdge),     yy(idx_TopEdge),     '  Top Edge',     'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
text(xx(idx_AboveCenter), yy(idx_AboveCenter), '  Above Center', 'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
text(xx(idx_Center),      yy(idx_Center),      '  Center',       'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
text(xx(idx_BelowCenter), yy(idx_BelowCenter), '  Below Center', 'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
text(xx(idx_BottomEdge),  yy(idx_BottomEdge),  '  Bottom Edge',  'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');
text(xx(idx_Node1),       yy(idx_Node1),       '  Node 1',       'Color', 'r', 'FontSize', 9, 'FontWeight', 'bold');

xlabel('x [m]');
ylabel('y [m]');
title('Tracked Nodes Locations (Red) vs All Nodes (Blue)');
axis equal;
grid on;
hold off;

% ------------------------------------------------------------------------
%                           Final outputs
% ------------------------------------------------------------------------

if flag_FinalPlots == 0

elseif flag_FinalPlots == 1
    
    % Find number of plots
    [s1,~] = size(PlotSettings);

    % Loop over plots
    for nplot = 1:s1

        % Figure number
        nfig = nplot;

        % Read plot field name
        cnodes_name = PlotSettings{nplot,1};

        % Check if compute damage for plotting
        if flag_BB == 0

        elseif flag_BB == 1
            if strcmp(cnodes_name,'Damage')
                % Compute damage
                phi = 1 - sum((u_NA>0).*V_NA,2)./phiD;
            end
        else
            error('flag_BB should be 0 or 1.')
        end

        % Read plot field variable
        cnodes = eval(PlotSettings{nplot,2});

        % Read plot settings
        ctitle        = PlotSettings{nplot,3};       % Colorbar title
        psize         = PlotSettings{nplot,4};       % Point size
        climits       = PlotSettings{nplot,5};       % Colormap limits
        cmap          = PlotSettings{nplot,6};       % Colormap
        box           = PlotSettings{nplot,7};       % Axes limits
        configuration = PlotSettings{nplot,8};       % Configuration: 'Reference' or 'Current'

        % Plot field
        if strcmp(configuration,'Reference')
            PlotField(nfig,xx,yy,cnodes,ctitle,psize,climits,cmap,box)
        elseif strcmp(configuration,'Current')
            PlotField(nfig,xx+v,yy+w,cnodes,ctitle,psize,climits,cmap,box)
        else
            error('Invalid configuration.')
        end

        % Save figure
        filename = ['../Outputs/' InputDeck '/' cnodes_name '.eps'];
        saveas(gcf,filename,'epsc');

    end

else

    error('flag_FinalPlots should be 0 or 1.')

end
