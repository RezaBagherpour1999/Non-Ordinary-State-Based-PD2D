% =========================================================================
% NOSB-PD2D: Non-Ordinary State-Based Peridynamics in MATLAB
%
% Original work Copyright (c) 2022, Oak Ridge National Laboratory (ORNL)
% Modified work Copyright (c) 2026 by Reza Bagherpour
%
% Description:
% NOSB-PD Force Density formulation, including deformation gradient,
% stress calculation, and zero-energy mode control.
%
% SPDX-License-Identifier: BSD-3-Clause
% =========================================================================

% ========================================================================
% The function ForceEnergyDensity computes the internal force density, 
% the strain energy density, and the bond stress metric for all nodes in 
% a 2D Non-Ordinary State-Based Peridynamics (NOSB-PD) formulation.
% ========================================================================

% Input
% -----
% xx                    : x-coordinates of all nodes in the reference configuration [Nnodes x 1]
% yy                    : y-coordinates of all nodes in the reference configuration [Nnodes x 1]
% v                     : displacement of each node in the x-direction [Nnodes x 1]
% w                     : displacement of each node in the y-direction [Nnodes x 1]
% E                     : Young's modulus
% nu                    : Poisson's ratio
% u_NA                  : neighbor node indices matrix [Nnodes x zmax]
% del                   : peridynamic horizon size (delta)
% IF_NA                 : influence function / bond status matrix [Nnodes x zmax]
%                         (1: intact bond, 0: broken bond)
% V_NA                  : reference area/volume of neighbor nodes [Nnodes x zmax]
% x_hat_NA              : initial x-coordinates of neighbor nodes/quadrature points [Nnodes x zmax]
% y_hat_NA              : initial y-coordinates of neighbor nodes/quadrature points [Nnodes x zmax]
% flag_RDUG             : if == 1, anti-symmetric pairwise force assembly is used
%                         for computational efficiency (Newton's 3rd law on uniform grids)
%                         if == 0, full loop assembly over all neighbors is used

% Output
% ------
% Fv                    : x-component of internal force density vector for all nodes [Nnodes x 1]
% Fw                    : y-component of internal force density vector for all nodes [Nnodes x 1]
% W                     : strain energy density for all nodes [Nnodes x 1]
% sigma_bond_stretch_ui : equivalent stress matrix (averaged maximum principal stress)
%                         evaluated for neighbor bonds [Nnodes x zmax]

% Discussion:
% ----------
% This function implements a 2D plane-stress Non-Ordinary State-Based Peridynamic 
% (NOSB-PD) correspondence model. The deformation gradient (F) is calculated via 
% the non-local shape tensor (K), from which the Green-Lagrange strain tensor (E), 
% Second Piola-Kirchhoff stress (S), and First Piola-Kirchhoff stress (P) are obtained.
%
% Zero-energy mode (hourglass) instability is controlled using an additional 
% stabilization state based on the formulation presented in:
%
% P. Li, Z. M. Hao, and W. Q. Zhen, "A stabilized non-ordinary state-based 
% peridynamic model," Computer Methods in Applied Mechanics and Engineering, 
% 339 (2018): 262–280.
%
% Reference for NOSB-PD correspondence framework:
% S. A. Silling, M. Epton, O. Weckner, J. Xu, and E. Askari, "Peridynamic states 
% and constitutive modeling," Journal of Elasticity, 88(2) (2007): 151–184.


function [Fv, Fw, W, sigma_bond_stretch_ui] = ForceEnergyDensity(xx, yy, v, w, E, nu, u_NA, del, IF_NA, V_NA, x_hat_NA, y_hat_NA,flag_RDUG)


% FORCEENERGYDENSITY Internal force and energy density for 2D NOSB-PD.
% 
% Description:
%   Calculates internal force density (Fv, Fw), macroelastic strain energy 
%   density (W), and equivalent bond stress for 2D Non-Ordinary State-Based 
%   Peridynamics (NOSB-PD). 
%
% Features:
%   - Fully vectorized kinematics and stress tensor computations for high performance.
%   - Analytic 2x2 tensor inversion with singularity regularization.
%   - Zero-energy mode instability control (Hourglass control) based on:
%     P. Li, Z.M. Hao, W.Q. Zhen, "A stabilized non-ordinary state-based 
%     peridynamic model", 2018.
%
% Inputs:
%   xx, yy        : Nodal coordinates (initial configuration).
%   v, w          : Nodal displacements in x and y directions.
%   c             : Peridynamic constant.
%   E, nu         : Young's Modulus and Poisson's ratio.
%   u_NA          : Neighbor index matrix (Nnodes x MaxNeighbors).
%   IF_NA         : Influence function / Bond status (1: intact, 0: broken).
%   V_NA          : Volume of neighboring nodes.
%   r_hat_NA      : Initial bond lengths.
%   x_hat_NA, y_hat_NA : Initial coordinates of neighboring nodes.
%
% Outputs:
%   Fv, Fw        : Internal force density vectors in x and y directions.
%   W             : Strain energy density for each node.
%   sigma_bond_stretch_ui : Equivalent stress matrix for bond failure criteria.

    % -------------------------- Setup & Initialization ------------------
    % NOTE: 'flag_RDUG' is used in the assembly loop but is missing from the 
    % original inputs. Uncomment the line below to define it locally and 
    % prevent the "Undefined variable" error without changing the signature.
    % flag_RDUG = 1; % 1: Reduced Data Usage (anti-symmetric), 0: Full loop

    Nnodes = length(xx);
    Fv = zeros(Nnodes, 1);
    Fw = zeros(Nnodes, 1);
    
    zmax = size(u_NA, 2);
    sigma_bond_stretch_ui = zeros(Nnodes, zmax);

    % Plane stress constitutive matrix components
    C11_mat = E / (1 - nu^2);
    C12_mat = nu * E / (1 - nu^2);
    C33_mat = E / (2 * (1 + nu)); % Shear Modulus (mu)

    % Pre-allocate vectorized components for 2x2 tensors (flattened to column vectors)
    K11 = zeros(Nnodes, 1); K12 = zeros(Nnodes, 1); 
    K21 = zeros(Nnodes, 1); K22 = zeros(Nnodes, 1);
    U11 = zeros(Nnodes, 1); U12 = zeros(Nnodes, 1); 
    U21 = zeros(Nnodes, 1); U22 = zeros(Nnodes, 1);

    % ====================================================================
    % STAGE 1: Vectorized Shape Tensor (K) & Deformation State (Uloc)
    % ====================================================================
    % Vectorizing over the neighbor index 'z' eliminates inner node loops 
    % and significantly reduces computational overhead.
    for z = 1:zmax
        uk = u_NA(:, z);
        valid = (uk > 0); % Process only existing neighbors
        
        if any(valid)
            ui_v = find(valid);
            uk_v = uk(valid);
            
            % Initial relative position vectors
            Dx = x_hat_NA(ui_v, z) - xx(ui_v);
            Dy = y_hat_NA(ui_v, z) - yy(ui_v);
            
            % Relative displacement vectors
            Dv = v(uk_v) - v(ui_v);
            Dw = w(uk_v) - w(ui_v);
            
            % Weighting factor (Influence Function * Node Volume)
            w_V = IF_NA(ui_v, z) .* V_NA(ui_v, z);
            
            % Accumulate Shape Tensor (K) components
            K11(ui_v) = K11(ui_v) + w_V .* Dx .* Dx;
            K12(ui_v) = K12(ui_v) + w_V .* Dx .* Dy;
            K21(ui_v) = K21(ui_v) + w_V .* Dy .* Dx;
            K22(ui_v) = K22(ui_v) + w_V .* Dy .* Dy;
            
            % Accumulate Relative Displacement Tensor components
            U11(ui_v) = U11(ui_v) + w_V .* Dv .* Dx;
            U12(ui_v) = U12(ui_v) + w_V .* Dv .* Dy;
            U21(ui_v) = U21(ui_v) + w_V .* Dw .* Dx;
            U22(ui_v) = U22(ui_v) + w_V .* Dw .* Dy;
        end
    end
    
    % ====================================================================
    % STAGE 2: Vectorized F, E, S, P, and K_inv (No for-loops)
    % ====================================================================
    
    % 1. Analytic regularization and inverse of Shape Tensor (K)
    eps_reg = 1e-14; % Small regularization term to prevent singularities
    K11 = K11 + eps_reg;
    K22 = K22 + eps_reg;
    
    detK = K11.*K22 - K12.*K21;
    invDetK = 1.0 ./ detK;
    
    Kinv11 =  K22 .* invDetK;
    Kinv22 =  K11 .* invDetK;
    Kinv12 = -K12 .* invDetK;
    Kinv21 = -K21 .* invDetK;
    
    % 2. Deformation Gradient Tensor (F = Uloc * Kinv + I)
    F11 = U11.*Kinv11 + U12.*Kinv21 + 1.0;
    F12 = U11.*Kinv12 + U12.*Kinv22;
    F21 = U21.*Kinv11 + U22.*Kinv21;
    F22 = U21.*Kinv12 + U22.*Kinv22 + 1.0;
    
    % 3. Green-Lagrange Strain Tensor (E = 0.5 * (F^T * F - I))
    E11 = 0.5 * (F11.*F11 + F21.*F21 - 1.0);
    E22 = 0.5 * (F12.*F12 + F22.*F22 - 1.0);
    E12 = 0.5 * (F11.*F12 + F21.*F22);
    
    % 4. Second Piola-Kirchhoff Stress Tensor (S = C:E) for Plane Stress
    S11 = C11_mat .* E11 + C12_mat .* E22;
    S22 = C12_mat .* E11 + C11_mat .* E22;
    S12 = 2 * C33_mat .* E12;
    
    % 5. First Piola-Kirchhoff Stress Tensor (P = F * S)
    P11 = F11.*S11 + F12.*S12;
    P12 = F11.*S12 + F12.*S22;
    P21 = F21.*S11 + F22.*S12;
    P22 = F21.*S12 + F22.*S22;
    
    % 6. Strain Energy Density (W = 0.5 * S:E)
    W = 0.5 * (S11.*E11 + S22.*E22 + 2 .* S12.*E12);
    
    % 7. Maximum Principal Stress Computation (p1)
    detF = F11.*F22 - F12.*F21;
    invDetF = 1.0 ./ detF;
    
    % Cauchy Stress Tensor components
    sig11 = (P11.*F11 + P12.*F12) .* invDetF;
    sig22 = (P21.*F21 + P22.*F22) .* invDetF;
    sig12 = (P11.*F21 + P12.*F22) .* invDetF;
    
    % First (maximum) principal stress formula
    p1 = 0.5.*(sig11+sig22) + sqrt(0.25.*(sig11-sig22).^2 + sig12.^2);
    
    % 8. Precompute P * Kinv (Critical for optimization in the force loop)
    PK11 = P11.*Kinv11 + P12.*Kinv21;
    PK12 = P11.*Kinv12 + P12.*Kinv22;
    PK21 = P21.*Kinv11 + P22.*Kinv21;
    PK22 = P21.*Kinv12 + P22.*Kinv22;
    
    % ====================================================================
    % STAGE 3: Pairwise Force Accumulation with Stabilization
    % ====================================================================
    
    % Penalty parameter for zero-energy mode control (Hourglass control)
    % Note: Horizon size (del) and constants are derived based on the specific grid.
    c_hg = 0.5 * (9 * E) / (pi * (del)^3); 
    
    for ui = 1:Nnodes
        xi = xx(ui); yi = yy(ui);
        vi = v(ui);  wi = w(ui);
        
        % Extract local tensors for node i
        F11_i = F11(ui); F12_i = F12(ui); F21_i = F21(ui); F22_i = F22(ui);
        PK11_i = PK11(ui); PK12_i = PK12(ui); PK21_i = PK21(ui); PK22_i = PK22(ui);
        
        for z = 1:zmax
            uk = u_NA(ui, z);
            
            % Assembly logic based on flag_RDUG (anti-symmetric optimization)
            % NOTE: 'flag_RDUG' must be defined (e.g., flag_RDUG = 1) or passed 
            % from the main script to avoid undefined variable errors.
            if (flag_RDUG == 1 && uk > ui) || (flag_RDUG == 0 && uk > 0)
                IFk = IF_NA(ui, z);
                if IFk == 0, continue; end % Skip broken bonds
                
                VkA = V_NA(ui, z);
                
                % Initial and deformed relative position vectors
                Dx = x_hat_NA(ui, z) - xi;
                Dy = y_hat_NA(ui, z) - yi;
                DX = (v(uk) - vi) + Dx;
                DY = (w(uk) - wi) + Dy;
                
                % Stabilization coefficient
                r2 = Dx*Dx + Dy*Dy;
                coef = c_hg / (r2 * sqrt(r2));
                
                C11_c = Dx*Dx; C12_c = Dx*Dy; C22_c = Dy*Dy;
                
                % Force state evaluation from node i's perspective
                ZI1 = DX - (F11_i*Dx + F12_i*Dy);
                ZI2 = DY - (F21_i*Dx + F22_i*Dy);
                
                TI1 = IFk * ((PK11_i*Dx + PK12_i*Dy) + coef * (C11_c*ZI1 + C12_c*ZI2));
                TI2 = IFk * ((PK21_i*Dx + PK22_i*Dy) + coef * (C12_c*ZI1 + C22_c*ZI2));
                
                % Force state evaluation from node k's perspective
                F11_k = F11(uk); F12_k = F12(uk); F21_k = F21(uk); F22_k = F22(uk);
                PK11_k = PK11(uk); PK12_k = PK12(uk); PK21_k = PK21(uk); PK22_k = PK22(uk);
                
                ZJ1 = -DX + F11_k*Dx + F12_k*Dy;
                ZJ2 = -DY + F21_k*Dx + F22_k*Dy;
                
                TJ1 = IFk * (-(PK11_k*Dx + PK12_k*Dy) + coef * (C11_c*ZJ1 + C12_c*ZJ2));
                TJ2 = IFk * (-(PK21_k*Dx + PK22_k*Dy) + coef * (C12_c*ZJ1 + C22_c*ZJ2));
                
                % Internal force density contribution
                dT1 = (TI1 - TJ1) * VkA;
                dT2 = (TI2 - TJ2) * VkA;
                
                % Accumulate forces for node i
                Fv(ui) = Fv(ui) + dT1;
                Fw(ui) = Fw(ui) + dT2;
                
                % Anti-symmetric reciprocal addition for node k (Newton's 3rd law)
                if flag_RDUG == 1
                    Fv(uk) = Fv(uk) - dT1;
                    Fw(uk) = Fw(uk) - dT2;
                end
            end
        end
    end
    
    % ====================================================================
    % STAGE 4: Vectorized Bond Stretch Evaluator (Principal Stress Based)
    % ====================================================================
    % Computes the equivalent stress metric used for the bond failure criteria 
    % by averaging the maximum principal stress of the two connected nodes.
    for z = 1:zmax
        uk = u_NA(:, z);
        valid = (uk > 0);
        
        if any(valid)
            ui_v = find(valid);
            uk_v = uk(valid);
            
            % p1 is the maximum principal stress; assigned directly for failure check
            sigma_bond_stretch_ui(ui_v, z) = 0.5 * (p1(ui_v) + p1(uk_v));
        end
    end

end
