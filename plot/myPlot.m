function myPlot(style, ns, qsz_all_pp,qsz_all_pp_label, qsz_all_pp_det,qsz_all_pp_det_label, qsz_all_dp, yax_label, ax_range)
    if strcmp(style, 'lineplot')
        myLinePlot(ns, qsz_all_pp,qsz_all_pp_label, qsz_all_pp_det,qsz_all_pp_det_label, qsz_all_dp, yax_label, ax_range);
    elseif strcmp(style, 'ratio')
        myRatioHeatmap(ns, qsz_all_pp, qsz_all_dp);
    else
    end
end

function myLinePlot(ns, qsz_all_pp,qsz_all_pp_label, qsz_all_pp_det,qsz_all_pp_det_label, qsz_all_dp, yax_label, ax_range)
    fontsize = 20;
    lindwidth = 2;
    markersize = 5;
    legend_pos = 'northwest';
    markers = ['o','+','*','s','x'];
    hold on; grid on;
    plot(ns.rho_all, qsz_all_pp,'r','linewidth',lindwidth,'DisplayName',strcat(' ',qsz_all_pp_label));
    plot(ns.rho_all, qsz_all_pp_det,'b','linewidth',lindwidth,'DisplayName',strcat(' ',qsz_all_pp_det_label));
    for i_eps = 1:ns.n_eps
        eps = ns.eps_all(i_eps);
        plot(ns.rho_all, qsz_all_dp(i_eps,:),...
            'linestyle','--','linewidth',lindwidth,...
            'marker',markers(i_eps),'markersize',markersize,...
            'DisplayName',sprintf(' DPS (\\epsilon=%s)',num2str(eps)));
    end
    legend('Location',legend_pos);
    set(gca, 'yscale', 'log', 'FontSize',fontsize);
    xlabel('$\mathbf{\rho}$','FontSize',fontsize+5,'FontWeight','bold','Interpreter','latex');
    if strcmp(yax_label,'qsz')
        ylabel('$\mathbf{\bar{Q}_t}$','FontSize',fontsize,'FontWeight','bold','Interpreter','latex');
%         ylabel('$\mathbf{E[Q_t]}$','FontSize',fontsize,'FontWeight','bold','Interpreter','latex');
    elseif strcmp(yax_label,'delay')
        ylabel('$\mathbf{\bar{W}_p}$','FontSize',fontsize,'FontWeight','bold','Interpreter','latex');
    end
    set(get(gca,'ylabel'), 'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right');
    
    cv = round(ns.coeff_var,2);
    ibr = round(ns.in_byte_rate,2);
    ar = round(ns.arrival_rate,2);
    title(strcat('$\mathbf{B_{in}}$=',num2str(ibr),', $\mathbf{\Lambda}$=',num2str(ar)),...
         'FontSize',fontsize+5,'FontWeight','bold','Interpreter','latex');

%     title(strcat('$\mathbf{s}$=',num2str(ns.scale), ', $\mathbf{B_{in}}$=',num2str(ibr),...
%         ', $\mathbf{\Lambda}$=',num2str(ar)),'FontSize',fontsize+5,'FontWeight','bold','Interpreter','latex');
%     title(strcat(ns.input_typ,', $B_{in}$=',num2str(ibr),', CV=',num2str(cv),...
%         ', $\Lambda$=',num2str(ar)),'Interpreter','latex');
    axis(ax_range);
    hold off;
end

function myRatioHeatmap(ns, qsz_all_pp, qsz_all_dp)
    try
        Z = qsz_all_pp'*ones(1,ns.n_eps);
    catch
        Z = qsz_all_pp*ones(1,ns.n_eps);
    end
    ratio = fliplr(qsz_all_dp./Z');
    h = heatmap(fliplr(round(ns.rho_all,2)), ns.eps_all, round(ratio,2)); % round(-log10(ratio),2)
    fmt = 'A = [%d';
    for n = 1:ns.N-1
        fmt = strcat(fmt, ', %d');
    end
    fmt = strcat(fmt, ']');
    h.Title = sprintf(fmt, ns.in_sizes);
%     strcat('-log_{10}(Q^{','DPS','}/Q^{',qsz_all_pp_label,'})')
    h.XLabel = '\rho';
    h.YLabel = '\epsilon';
%     h.ColorLimits = [-1 5];
    h.ColorLimits = [0 1];
    h.FontSize = 20;
    colormap summer;
end