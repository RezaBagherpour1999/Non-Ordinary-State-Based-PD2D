
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
% Input deck for a crack branching problem on a pre-notched soda–lime glass 
% thin plate 
% ========================================================================

% Reference: 
% ---------
% F. Bobaru and G. Zhang, Why do cracks branch? A peridynamic investigation of
% dynamic brittle fracture, International Journal of Fracture 196 (2015): 59-98.

% ------------------------------------------------------------------------
% 1. Domain Geometry and Spatial Discretization
% ------------------------------------------------------------------------

% Domain bounding coordinates [m]
% Plate dimensions: Length L = 0.10 m (100 mm), Height H = 0.04 m (40 mm)
Xo = -0.05; % [m] : Minimum x-coordinate (left edge)
Xn =  0.05; % [m] : Maximum x-coordinate (right edge)
Yo = -0.02; % [m] : Minimum y-coordinate (bottom edge)
Yn =  0.02; % [m] : Maximum y-coordinate (top edge)

% Grid resolution (number of nodal divisions)
Nx = 300;   % Number of nodes along the x-direction
Ny = 120;   % Number of nodes along the y-direction

% Spatial nodal spacing along x-axis [m] (dx = 0.10 / 300 = 0.333 mm)
dx = (Xn - Xo) / Nx;

% Grid perturbation parameter:
% PG = 0 enforces an ideal uniform lattice; PG > 0 introduces random nodal jitter
PG = 0;

% ------------------------------------------------------------------------
% 2. Temporal Discretization & Time-Integration Solver
% ------------------------------------------------------------------------

% Temporal range [s]
Ti = 0.0;     % Initial simulation time [s]
Tf = 4.0e-5;  % Final simulation time [s] (40.0 microseconds)

% Stable time step size [s] (must satisfy CFL condition: dt < dx / c_wave)
dt = 5.0e-8;  % Discrete time increment [s] (20.0 nanoseconds = 0.02 microseconds)

% Numerical explicit time-stepping scheme
TimeScheme = 'VVerlet'; % Velocity-Verlet algorithm

% Artificial viscous damping coefficient.
% Set to 0 to disable artificial damping (for pure elastodynamics/dynamic fracture).
% Note: For quasi-static problems, a high damping value (e.g., 8e8) should be 
% applied, scaled appropriately according to the total simulation time.
c_damp = 0;

% ------------------------------------------------------------------------
% 3. Peridynamic (PD) Constitutive Model
% ------------------------------------------------------------------------

% Peridynamic formulation type
model = 'NOSB-PD'; % Generalized Prototype Microelastic Brittle (GPMB) model

% Kinematic state assumption for 2D analysis
PlanarModel = 'PlaneStress'; % Plane stress formulation for thin plates

% Non-local horizon radius [m] (delta = m * dx, with m = 3.0)
del = dx * 3.0;

% Radial influence function exponent (omega = 0 corresponds to constant kernel)
omega = 0;

% Bond failure control switch:
% flag_BB = 1 enables irreversible bond breaking; 0 represents purely elastic response
flag_BB = 1;

% ------------------------------------------------------------------------
% 4. Material Properties (Soda–Lime Glass)
% ------------------------------------------------------------------------

% Mass density [kg/m^3]
rho = 2235;

% Elastic mechanical properties
% Note: Soda-lime glass Young's modulus typically spans 65–72 GPa
E  = 65.0e+9; % Modulus of elasticity [Pa] (65 GPa)
nu = 0.20;    % Poisson's ratio [-]

% Critical fracture energy / critical energy release rate [J/m^2]
Go = 3.8;     % Fracture energy [J/m^2]

% ------------------------------------------------------------------------
% 5. Meshfree Horizon Integration Method
% ------------------------------------------------------------------------

% Volume/area integration correction algorithm for horizon intersection
% 'FA': Fractional Area method (exact geometric boundary intersection)
AlgName = 'FA';

% ------------------------------------------------------------------------
% 6. Boundary Conditions & External Loading
% ------------------------------------------------------------------------

% Spatial nodal spacing along y-axis [m] (dy = 0.04 / 120 = 0.333 mm)
dy = (Yn - Yo) / Ny;

% External traction magnitude applied at top and bottom boundaries [Pa]
sigma = 2.0e6; % 2.0 MPa tensile stress

% Equivalent body force density functions [N/m^3]:
% Peridynamic boundary traction is converted into an equivalent fictitious body 
% force density spread across the boundary node layer of thickness dy: b = sigma / dy.
bvfunc = @(x,y,t) (0.*x + 0.*y) * t;                          % x-component: Zero horizontal body force
bwfunc = @(x,y,t) ( abs(y) > Yn - dy ) .* sigma .* sign(y)/dy; % y-component: Symmetric tensile force density

% ------------------------------------------------------------------------
% 7. Initial Kinematic Conditions
% ------------------------------------------------------------------------

% Initial nodal displacement fields [m] (quiescent initial state)
vofunc = @(x,y) (0.*x + 0.*y); % Initial displacement along x (horizontal)
wofunc = @(x,y) (0.*x + 0.*y); % Initial displacement along y (vertical)

% Initial nodal velocity fields [m/s]
Vvofunc = @(x,y) (0.*x + 0.*y); % Initial velocity along x
Vwofunc = @(x,y) (0.*x + 0.*y); % Initial velocity along y

% ------------------------------------------------------------------------
% 8. No-Fail Boundary Layer Zone
% ------------------------------------------------------------------------

% Nodal mask function to suppress non-physical bond breaking near external load boundaries.
% Bonds connected to nodes within a boundary band of width equal to horizon (delta) are preserved.
nofailfunc = @(x,y) ( abs(y) > Yn - del );

% ------------------------------------------------------------------------
% 9. Initial Crack (Pre-Notch) Geometry
% ------------------------------------------------------------------------

% Coordinates defining the initial notch line segment: [Xc1, Yc1, Xc2, Yc2]
% Horizontal slit along the mid-plane (y = 0) from left edge (x = -0.05) to domain center (x = 0.0)
% Pre-crack length a0 = 0.05 m (50 mm)
% Note: For an uncracked domain (crack-free problem), define the pre-crack segment
% completely outside the computational domain (e.g., [-1.0, -1.0, -1.0, -1.0]).

PreNotchCoordinates = [-0.05, 0.0, 0.0, 0.0];

% ------------------------------------------------------------------------
% 10. Post-Processing & Runtime Visualization
% ------------------------------------------------------------------------

% Dynamic visualization during explicit time-stepping
flag_DynamicPlotting = 1;      % 1 = Enable real-time dynamic plots, 0 = Disable

% Output plot frequency during simulation
DynamicPlotFrequency = 40;     % Render graphical frame every 40 time steps

% Command line solver progress display frequency
TimeStepDisplayFrequency = 20; % Print step status every 20 time increments

% Final state visualization
flag_FinalPlots = 1;           % 1 = Generate final solution plots at t = Tf

% Field visualization configuration table:
% Format: {Field Name, Field Variable, Colorbar Label, Point Size, Colormap Limits, Colormap, Axis Limits, Configuration}
PlotSettings = { ...
    'StrainEnergyDensity', 'log10(W)', '$\log_{10}(W)$', 8, [0 3.5], 'jet'   , [Xo Xn Yo Yn], 'Reference'; ...
    'Damage'             , 'phi'     , '$\varphi$'      , 8, [0 0.4], 'parula', [Xo Xn Yo Yn], 'Reference'  ...
};

% Visualization flag for pre-notch representation
% 1 = Display initial pre-notch surface as fully damaged bonds (phi = 1)
flag_DamagedPrenotches = 1;

% ========================================================================
% End of Input Deck
% ========================================================================
---------------------------------------
