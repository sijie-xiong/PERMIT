
% shaper = 'DPS_VSCI';
shaper = 'VSCI_CSCI';

subfolder = 'IoT';
rp = ns.get_rp(subfolder);
res_prefix = strcat(shaper, '_opt');

ns = NetSetting('where','local','synthetic',false,...
                'device','Sleep','shaper',shaper);

eps_s_actual_all = NaN(ns.n_eps, ns.n_rho);
eps_t_actual_all = NaN(ns.n_eps, ns.n_rho);
rho_est_all = NaN(ns.n_eps, ns.n_rho);
qsz_est_all = NaN(ns.n_eps, ns.n_rho);
Q_opt_all = NaN(ns.N, ns.hN, ns.n_eps, ns.n_rho);

for i_eps = 1:ns.n_eps
    for i_rho = 1:ns.n_rho
        ns.i_eps = i_eps-1; ns.i_rho = i_rho-1;
        fn = ns.get_fn(res_prefix, 'ind');
        try
            opt_res = load(strcat(rp, '/', fn));
        catch
            continue;
        end
        eps_s_actual_all(i_eps, i_rho) = opt_res.eps_s_actual;
        eps_t_actual_all(i_eps, i_rho) = opt_res.eps_t_actual;
        rho_est_all(i_eps, i_rho) = opt_res.rho_est;
        qsz_est_all(i_eps, i_rho) = opt_res.qsz_est;
        Q_opt_all(:,:,i_eps, i_rho) = opt_res.Q_opt;
    end
end

res = struct();
res.ns = ns;
res.eps_s_actual_all = eps_s_actual_all;
res.eps_t_actual_all = eps_t_actual_all;
res.rho_est_all = rho_est_all;
res.qsz_est_all = qsz_est_all;
res.Q_opt_all = Q_opt_all;

% Local (save to a local folder)
fn = ns.get_fn(res_prefix, 'agg');
save(strcat(rp,'/',fn), '-struct', 'res');