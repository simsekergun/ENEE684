function [nodes, elems, epsilon_r] = build_simple_mesh(varargin)
% BUILD_SIMPLE_MESH  Conforming triangular mesh for a rectangular dielectric
%                    core embedded in a rectangular cladding domain.
%
%   The four edges of the waveguide core rectangle are inserted as
%   CONSTRAINED edges so that no triangle ever straddles the core/cladding
%   boundary.  Every triangle is purely inside or purely outside the core.
%
%   KEY FIX vs previous version:
%     Constraint indices are recorded AT THE TIME points are added to the
%     list, so they are exact integers — no snap-tolerance searching needed.
%     uniquetol is NOT called on the full point set because it can silently
%     shift interface points and break constraint chains.  Instead, corners
%     (which appear on two segments) are deduplicated explicitly.
%
% ── Syntax ────────────────────────────────────────────────────────────────
%   [nodes, elems, epsilon_r] = build_simple_mesh()
%   [nodes, elems, epsilon_r] = build_simple_mesh(Name, Value, ...)
%
% ── Parameters ────────────────────────────────────────────────────────────
%   'w_core'      core width          (µm)  default 0.890
%   'h_core'      core height         (µm)  default 0.670
%   'w_sim'       simulation width    (µm)  default 6.0
%   'h_sim'       simulation height   (µm)  default 4.0
%   'n_core'      core refractive index     default 1.9955  (Si3N4 @1550 nm)
%   'n_clad'      clad refractive index     default 1.4440  (SiO2  @1550 nm)
%   'target_edge_length'  target edge length  (µm)  default 0.15
%   'plot_mesh'   show mesh figure          default false
%
% ── Outputs ───────────────────────────────────────────────────────────────
%   nodes      [Nn x 2]  node (x,y) coordinates in µm
%   elems      [Ne x 3]  element connectivity (1-based)
%   epsilon_r  [Ne x 1]  relative permittivity per element

%% ── Parse inputs ──────────────────────────────────────────────────────────
p = inputParser;
addParameter(p, 'w_core',     0.890,  @isnumeric);
addParameter(p, 'h_core',     0.670,  @isnumeric);
addParameter(p, 'w_sim',      6.0,    @isnumeric);
addParameter(p, 'h_sim',      4.0,    @isnumeric);
addParameter(p, 'n_core',     1.9955, @isnumeric);
addParameter(p, 'n_clad',     1.4440, @isnumeric);
addParameter(p, 'target_edge_length', 0.15,   @isnumeric);
addParameter(p, 'plot_mesh',  false,  @islogical);
parse(p, varargin{:});
o = p.Results;

edge_length  = o.target_edge_length;
min_dist = edge_length * 0.25;
gap = edge_length * 0.45;
eps_box = edge_length * 1e-4;   % tiny tolerance only for the outer box keep/discard
eps_core = edge_length * 1e-4;
tol_v = edge_length * 1e-4;

% Outer box
xL = -o.w_sim/2;   xR =  o.w_sim/2;
yB = -o.h_sim/2;   yT =  o.h_sim/2;

% Core rectangle (horizontally centred)
xCL = -o.w_core/2;  xCR =  o.w_core/2;
yCB =  -o.h_core/2;          yCT =  o.h_core/2;

fprintf('[build_simple_mesh]  core %.4f x %.4f µm  |  box %.3f x %.3f µm  |  edge_length %.4f µm\n', ...
    o.w_core, o.h_core, o.w_sim, o.h_sim, edge_length);

%% ── Step 1: Build point list and constraint list TOGETHER ─────────────────
%   Strategy: Each constrained segment is sampled with linspace.  The indices of
%   those points in the global list are known exactly at insertion time,
%   so we can build constraint pairs [i, i+1] directly — no searching,
%   no snap tolerance, no way for the chain to break.
%
%   Corner sharing: Adjacent segments share one endpoint.  We handle this by always
%   accepting the start-point of the next segment to be the same index as
%   the end-point of the previous segment (no duplicate point inserted).

pts = zeros(0, 2);          % point list (grows below)
C   = zeros(0, 2, 'int32'); % constraint pairs

% Number of intervals on each edge (at least 2 so endpoints are included)
nOB = max(2, round(o.w_sim  / edge_length*o.n_clad/o.n_core));   % outer bottom/top
nOH = max(2, round(o.h_sim  / edge_length*o.n_clad/o.n_core));   % outer left/right
nCW = max(2, round(o.w_core / edge_length));   % core bottom/top
nCH = max(2, round(o.h_core / edge_length));   % core left/right

% ── helper: append a segment, returning updated pts and C ────────────────
%   x0,y0 -> x1,y1   with n+1 points (n intervals).
%   If reuse_start = true, the first point already exists at row idx0
%   and must NOT be re-added; idx0 is then the index to reuse.
    function [pts, C, idx_end] = append_seg(pts, C, x0,y0, x1,y1, n, reuse_start, idx0)
        t     = linspace(0, 1, n+1)';
        new_x = x0 + t*(x1-x0);
        new_y = y0 + t*(y1-y0);
        if reuse_start
            % first point already in pts at row idx0; skip it
            new_pts = [new_x(2:end), new_y(2:end)];
            first   = idx0;
        else
            new_pts = [new_x, new_y];
            first   = size(pts,1) + 1;
        end
        pts = [pts; new_pts];
        % indices of the n+1 points
        if reuse_start
            idxs = [idx0; (size(pts,1) - n + 1 : size(pts,1))'];
        else
            idxs = (first : first + n)';
        end
        % chain consecutive pairs
        pairs = int32([idxs(1:end-1), idxs(2:end)]);
        C     = [C; pairs];
        idx_end = idxs(end);
    end

% ── Outer box (closed loop: BL -> BR -> TR -> TL -> BL) ──────────────────
[pts, C, i_BR] = append_seg(pts, C,  xL,yB,  xR,yB,  nOB, false, []);
[pts, C, i_TR] = append_seg(pts, C,  xR,yB,  xR,yT,  nOH, true,  i_BR);
[pts, C, i_TL] = append_seg(pts, C,  xR,yT,  xL,yT,  nOB, true,  i_TR);
i_BL           = 1;   % bottom-left was the very first point added
[pts, C,    ~] = append_seg(pts, C,  xL,yT,  xL,yB,  nOH, true,  i_TL);
% close outer loop: last point of left edge must coincide with i_BL
% (linspace guarantees pts(end) == (xL,yB) == pts(1), but they are
%  separate rows — add one closing constraint)
C = [C; int32([size(pts,1), i_BL])];

% ── Core rectangle (closed loop: BL -> BR -> TR -> TL -> BL) ─────────────
[pts, C, ic_BR] = append_seg(pts, C,  xCL,yCB,  xCR,yCB,  nCW, false, []);
[pts, C, ic_TR] = append_seg(pts, C,  xCR,yCB,  xCR,yCT,  nCH, true,  ic_BR);
[pts, C, ic_TL] = append_seg(pts, C,  xCR,yCT,  xCL,yCT,  nCW, true,  ic_TR);
ic_BL           = ic_BR - nCW;   % first point of the core bottom edge
[pts, C,     ~] = append_seg(pts, C,  xCL,yCT,  xCL,yCB,  nCH, true,  ic_TL);
C = [C; int32([size(pts,1), ic_BL])];

% ── Interior scatter (no constraints, just fill points) ───────────────────
%   Keep scatter points strictly away from all boundaries to avoid
%   accidentally duplicating constrained points.
[gx, gy] = meshgrid( ...
    linspace(xL  + gap, xR  - gap, max(4, round(o.w_sim  / edge_length))), ...
    linspace(yB  + gap, yT  - gap, max(4, round(o.h_sim  / edge_length))));
scatter_pts = [gx(:), gy(:)];

% Remove scatter points that are too close to any constrained boundary
% (avoids degeneracy / near-duplicate nodes in the triangulation)
on_outer = abs(scatter_pts(:,1)-xL)  < min_dist | ...
           abs(scatter_pts(:,1)-xR)  < min_dist | ...
           abs(scatter_pts(:,2)-yB)  < min_dist | ...
           abs(scatter_pts(:,2)-yT)  < min_dist;
on_core  = abs(scatter_pts(:,1)-xCL) < min_dist | ...
           abs(scatter_pts(:,1)-xCR) < min_dist | ...
           abs(scatter_pts(:,2)-yCB) < min_dist | ...
           abs(scatter_pts(:,2)-yCT) < min_dist;
scatter_pts = scatter_pts(~on_outer & ~on_core, :);

pts = [pts; scatter_pts];

% Remove duplicate constraint pairs (same unordered pair listed twice)
C = unique(sort(double(C), 2), 'rows');

fprintf('                     seed points: %d   constraint edges: %d\n', ...
    size(pts,1), size(C,1));

%% ── Step 2: Constrained Delaunay triangulation ───────────────────────────
fprintf('                     running constrained Delaunay ...\n');
DT = delaunayTriangulation(pts, C);

%% ── Step 3: Discard triangles outside the simulation box ─────────────────
conn_all  = DT.ConnectivityList;
nodes_all = DT.Points;

cx_all = (nodes_all(conn_all(:,1),1) + nodes_all(conn_all(:,2),1) + nodes_all(conn_all(:,3),1)) / 3;
cy_all = (nodes_all(conn_all(:,1),2) + nodes_all(conn_all(:,2),2) + nodes_all(conn_all(:,3),2)) / 3;

keep  = cx_all >= xL - eps_box & cx_all <= xR + eps_box & ...
        cy_all >= yB - eps_box & cy_all <= yT + eps_box;
elems_keep = conn_all(keep, :);

%% ── Step 4: Compact node list ────────────────────────────────────────────
used     = unique(elems_keep(:));
node_map = zeros(size(nodes_all,1), 1);
node_map(used) = 1:numel(used);
nodes = nodes_all(used, :);
elems = node_map(elems_keep);

Nn = size(nodes, 1);
Ne = size(elems, 1);

%% ── Step 5: Assign permittivity by centroid ──────────────────────────────
%   Because the constraints guarantee that no triangle straddles the
%   boundary, the centroid test is exact — a centroid is either strictly
%   inside or strictly outside the core rectangle.
%   We use a tiny positive tolerance (1e-4 * edge_length) purely to handle
%   numerical rounding in the centroid coordinates.

cx = (nodes(elems(:,1),1) + nodes(elems(:,2),1) + nodes(elems(:,3),1)) / 3;
cy = (nodes(elems(:,1),2) + nodes(elems(:,2),2) + nodes(elems(:,3),2)) / 3;

in_core  = cx > xCL - eps_core & cx < xCR + eps_core & ...
           cy > yCB - eps_core & cy < yCT + eps_core;

epsilon_r              = ones(Ne,1) * o.n_clad^2;
epsilon_r(in_core)     = o.n_core^2;

n_c  = sum( in_core);
n_cl = sum(~in_core);
fprintf('                     nodes: %d   elems: %d   (core: %d  clad: %d)\n', Nn, Ne, n_c, n_cl);

%% ── Step 6: Verify — zero triangles straddle the core boundary ───────────
%   A triangle straddles if it has at least one node strictly inside AND
%   at least one node strictly outside (nodes exactly on the boundary are
%   excluded from both sets to avoid false positives).

strictly_inside  = nodes(:,1) > xCL + tol_v & nodes(:,1) < xCR - tol_v & ...
                   nodes(:,2) > yCB + tol_v & nodes(:,2) < yCT - tol_v;
strictly_outside = nodes(:,1) < xCL - tol_v | nodes(:,1) > xCR + tol_v | ...
                   nodes(:,2) < yCB - tol_v | nodes(:,2) > yCT + tol_v;

n1 = elems(:,1); n2 = elems(:,2); n3 = elems(:,3);
has_in  = strictly_inside(n1)  | strictly_inside(n2)  | strictly_inside(n3);
has_out = strictly_outside(n1) | strictly_outside(n2) | strictly_outside(n3);
straddle = has_in & has_out;

if any(straddle)
    warning('build_simple_mesh: %d triangle(s) STILL straddle the boundary.', sum(straddle));
    % Print the first few offending centroids to help debug
    bad_idx = find(straddle, 5);
    fprintf('  First offending centroids (x, y):\n');
    for k = bad_idx'
        fprintf('    elem %d:  cx=%.6f  cy=%.6f\n', k, cx(k), cy(k));
    end
else
    fprintf('                     VERIFIED: 0 triangles straddle the core boundary.\n');
end

%% ── Step 7: Optional mesh plot ───────────────────────────────────────────
if o.plot_mesh
    figure('Name','build_simple_mesh','Color','w','Position',[80 80 1000 620]);

    % Full mesh coloured by material
    ax1 = subplot(1,2,1);
    patch('Parent',ax1, 'Faces',elems(~in_core,:), 'Vertices',nodes, ...
          'FaceColor',[0.75 0.85 1.0], 'EdgeColor',[0.55 0.65 0.85], ...
          'LineWidth',0.25, 'DisplayName','Cladding');
    hold(ax1,'on');
    patch('Parent',ax1, 'Faces',elems(in_core,:), 'Vertices',nodes, ...
          'FaceColor',[1.0 0.78 0.35], 'EdgeColor',[0.80 0.50 0.05], ...
          'LineWidth',0.25, 'DisplayName','Core');
    xcore = [xCL xCR xCR xCL xCL];
    ycore = [yCB yCB yCT yCT yCB];
    plot(ax1, xcore, ycore, 'r-', 'LineWidth',1.8, 'DisplayName','Core boundary');
    axis(ax1,'equal','tight');
    xlabel(ax1,'x  (µm)'); ylabel(ax1,'y  (µm)');
    title(ax1, sprintf('Full mesh — %d nodes, %d elems', Nn, Ne));
    legend(ax1,'Location','northeast','FontSize',8);

    % Zoom to top-right corner of core
    ax2 = subplot(1,2,2);
    zm  = edge_length * 4;
    patch('Parent',ax2, 'Faces',elems(~in_core,:), 'Vertices',nodes, ...
          'FaceColor',[0.75 0.85 1.0], 'EdgeColor',[0.40 0.55 0.80], 'LineWidth',0.5);
    hold(ax2,'on');
    patch('Parent',ax2, 'Faces',elems(in_core,:), 'Vertices',nodes, ...
          'FaceColor',[1.0 0.78 0.35], 'EdgeColor',[0.75 0.42 0.00], 'LineWidth',0.5);
    plot(ax2, xcore, ycore, 'r-', 'LineWidth',2);
    % mark straddling triangles in magenta (should be none)
    if any(straddle)
        patch('Parent',ax2,'Faces',elems(straddle,:),'Vertices',nodes, ...
              'FaceColor','m','EdgeColor','k','LineWidth',1.0,'DisplayName','STRADDLE');
    end
    axis(ax2,'equal');
    xlim(ax2,[xCR-zm, xCR+zm]);
    ylim(ax2,[yCT-zm, yCT+zm]);
    xlabel(ax2,'x  (µm)'); ylabel(ax2,'y  (µm)');
    title(ax2,'Top-right corner zoom');
    box(ax2,'on');

    sgtitle(sprintf('Si_3N_4 core (n=%.4f) in SiO_2 (n=%.4f)', o.n_core, o.n_clad), 'Interpreter','tex');
end

end % ── build_simple_mesh ─

