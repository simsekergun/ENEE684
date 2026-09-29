%% ============================================================
%  LOCAL FUNCTION: plot_H_mode
%  Plot |Hx|, |Hy|, |Hz| and optionally the full E+H panel
%  for one mode.
%% ============================================================
function plot_EH_mode(md, m_idx, save_fig)

    nodes = md.nodes;
    elems = md.elems;

    Hx = md.Hx_centroid;
    Hy = md.Hy_centroid;
    Hz = md.Hz_centroid;

    n_eff_str = sprintf('%.6f + %.3ei', real(md.n_eff), imag(md.n_eff));

    % ── Figure 3: Side-by-side E and H (6 panels) ─────────────────────────
    if isfield(md, 'Ex_centroid')
        Ex = md.Ex_centroid;
        Ey = md.Ey_centroid;
        Ez = md.Ez_centroid;

        fig3 = figure('Name', sprintf('Mode %d  E and H fields', m_idx), ...
                      'Color', 'w', 'Position', [160, 160, 1400, 680]);

        all_cmps = {abs(Ex), abs(Ey), abs(Ez), abs(Hx), abs(Hy), abs(Hz)};
        all_lbl  = {'|E_x|','|E_y|','|E_z|','|H_x|','|H_y|','|H_z|'};

        for k = 1:6
            ax = subplot(2,3,k);
            tripcolor_patch(ax, nodes, elems, all_cmps{k});
            %clim(ax, [0, 1]);
            colorbar(ax);
            title(ax, all_lbl{k}, 'Interpreter','tex', 'FontSize',11);
            xlabel(ax, 'x  (µm)'); ylabel(ax, 'y  (µm)');
            axis(ax, 'equal', 'tight');
        end

        % Row labels
        annotation(fig3,'textbox',[0.01,0.72,0.05,0.08], ...
            'String','E-field','FontWeight','bold','EdgeColor','none', ...
            'FontSize',11,'Rotation',90,'HorizontalAlignment','center');
        annotation(fig3,'textbox',[0.01,0.22,0.05,0.08], ...
            'String','H-field','FontWeight','bold','EdgeColor','none', ...
            'FontSize',11,'Rotation',90,'HorizontalAlignment','center');

        sgtitle(fig3, ...
            sprintf('Mode %d   E and H field magnitudes   n_{eff} = %s', ...
            m_idx, n_eff_str), 'Interpreter','tex', 'FontSize',11);

        if save_fig
            fname = sprintf('mode%d_EH_fields.png', m_idx);
            print(fig3, fname, '-dpng', '-r200');
            fprintf('  Saved %s\n', fname);
        end
    end

end 
% ── plot_EH_mode ──

