%% phase_to_hole_diameter.m
% Ku-band ABS perforated GRIN lens synthesis: Fermat phase map -> local
% n_eff -> Maxwell-Garnett fill fraction -> drilled hole diameter grid.
%
% IMPORTANT SIGN NOTE (read this before trusting any output):
% phi_error(x,y) = k0*[sqrt(F_H^2+x^2)-F_H] + k0*[sqrt(F_E^2+y^2)-F_E]
% is the *spherical-wave phase error* relative to the on-axis point (zero
% at the center, growing toward the edges). A perforated ABS slab is a
% delay-only medium (n_eff >= 1, it can only slow the wave down, never
% speed it up). Since the on-axis ray has the SHORTEST physical path from
% the horn phase center, it is the one that needs the MOST extra delay to
% be brought back into phase with the (already-delayed) edge rays.
% => the physical lens delay is the COMPLEMENT of phi_error, not
%    phi_error itself. Mapping phi_error directly to n_eff (increasing
%    outward) builds a lens that is thickest/densest at the EDGES and
%    most porous at the CENTER -> that DEFOCUSES the beam instead of
%    correcting it. This script applies the correct inversion (see psi
%    below) and wraps automatically into Fresnel zones if the required
%    delay ever exceeds 360 deg (relevant if you scale this up to the
%    larger Q2 aperture).

clear; clc;

%% ---- User inputs ----
f0      = 15e9;         % design (center) frequency, Hz
c       = 3e8;
lambda0 = c/f0;
k0      = 2*pi/lambda0;

% Original antenna: throat 23x10mm, aperture 120x80mm -- these are your
% own originally-verified phase-center distances (F_H, F_E) from your
% first Ansys phase-map plots, matched exactly to your 221.2deg/99.6deg data.
F_H = 145.28e-3;        % H-plane phase-center distance to aperture, m
F_E = 149.10e-3;        % E-plane phase-center distance to aperture, m
x_half = 60e-3;         % aperture half-width, H-plane, m (120 mm aperture)
y_half = 40e-3;         % aperture half-width, E-plane, m (80 mm aperture)

eps_ABS   = 2.6;        % stock ABS -
tan_d_ABS = 0.008;      % typical stock ABS tan-delta; MEASURE on your own printed coupons regardless
f_max_porosity = 0.65;  % practical circular-hole area fraction ceiling
                         % (structural/print-resolution limit; theoretical
                         % max for a square lattice of circles is pi/4=0.785)

pitch = 3.0e-3;         % unit-cell pitch, mm -> must satisfy pitch <= lambda_min/4
                         % lambda_min at 18 GHz = 16.7 mm -> pitch <= ~4.2 mm, 3 mm is safe

d_min_print = 0.4e-3;   % smallest reliably printable hole diameter (FDM nozzle-dependent)

%% ---- 1. Rebuild the Fermat phase-error map (same math as your plots) ----
x = -x_half:pitch:x_half;
y = -y_half:pitch:y_half;
[X,Y] = meshgrid(x,y);

phi_H = k0*(sqrt(F_H^2 + X.^2) - F_H);
phi_E = k0*(sqrt(F_E^2 + Y.^2) - F_E);
phi_error = phi_H + phi_E;                 % rad, 0 at center, grows outward

%% ---- 2. Convert to the physical LENS DELAY the ABS must supply ----
C = max(phi_error(:));                     % reference = max error over aperture
psi = mod(C - phi_error, 2*pi);            % dome-shaped: max at center, ->0 at edges
                                            % (auto-wraps into extra Fresnel zones
                                            %  if the aperture/required range grows)

%% ---- 3. psi -> required effective index (uniform-thickness slab) ----
eps_i = 1;                                  % air
n_ABS = sqrt(eps_ABS);

eps_min = mg_forward(f_max_porosity, eps_ABS, eps_i);
n_min   = sqrt(eps_min);
n_max   = n_ABS;

t_lens = lambda0/(n_max - n_min);           % thickness giving a full 360 deg span
fprintf('Lens thickness (uniform): %.2f mm (%.2f * lambda0)\n', t_lens*1e3, t_lens/lambda0);

n_eff = n_min + (n_max-n_min).*psi/(2*pi);  % linear map, psi=0->n_min, psi=2pi(=0)->n_max

%% ---- 4. Invert Maxwell-Garnett (2D, cylindrical air holes) for fill fraction ----
eps_eff = n_eff.^2;
f_fill  = mg_inverse(eps_eff, eps_ABS, eps_i);
f_fill  = min(max(f_fill, 0), f_max_porosity);   % clip to physical/practical range

%% ---- 5. Fill fraction -> hole diameter on the chosen square lattice ----
d = pitch*sqrt(4*f_fill/pi);
d(d < d_min_print) = 0;                    % below this, print as solid (no hole)
d = min(d, 0.9*pitch);                     % keep a manufacturable wall between holes

%% ---- 6. Validate: recompute achieved phase from the CLIPPED design ----
eps_eff_real = mg_forward(f_fill, eps_ABS, eps_i);
n_eff_real   = sqrt(eps_eff_real);
psi_real     = k0*(n_eff_real-1)*t_lens;           % achieved delay, unwrapped
phase_error_deg = mod(rad2deg(psi_real - psi)+180,360)-180;
fprintf('Max residual phase error after clipping: %.1f deg (RMS %.1f deg)\n', ...
    max(abs(phase_error_deg(:))), rms(phase_error_deg(:)));

%% ---- 7. Plots ----
figure;
subplot(1,2,1);
surf(X*1e3, Y*1e3, rad2deg(psi), 'EdgeColor','none'); view(2); colorbar; axis equal tight;
title('Required lens delay \psi(x,y) [deg] (correct sign: dome, max at center)');
xlabel('x - H-plane (mm)'); ylabel('y - E-plane (mm)');

subplot(1,2,2);
surf(X*1e3, Y*1e3, d*1e3, 'EdgeColor','none'); view(2); colorbar; axis equal tight;
title('Hole diameter map (mm)');
xlabel('x - H-plane (mm)'); ylabel('y - E-plane (mm)');

%% ---- 8. Export hole table for CAD / print-path generation ----
T = table(X(:)*1e3, Y(:)*1e3, d(:)*1e3, 'VariableNames', {'x_mm','y_mm','d_mm'});
writetable(T, 'lens_hole_table.csv');
fprintf('Wrote %d hole positions to lens_hole_table.csv\n', height(T));

%% ---- Local functions: 2D Maxwell-Garnett (cylindrical inclusions, transverse field) ----
function eps_eff = mg_forward(f, eps_h, eps_i)
    % eps_eff = eps_h*[(eps_i+eps_h)+f*(eps_i-eps_h)] / [(eps_i+eps_h)-f*(eps_i-eps_h)]
    eps_eff = eps_h.*((eps_i+eps_h)+f.*(eps_i-eps_h)) ./ ((eps_i+eps_h)-f.*(eps_i-eps_h));
end

function f = mg_inverse(eps_eff, eps_h, eps_i)
    % Closed-form inverse of mg_forward for f (air fill fraction)
    f = (eps_h+eps_i).*(eps_eff-eps_h) ./ ((eps_i-eps_h).*(eps_h+eps_eff));
end
