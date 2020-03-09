
setenv('i_eps', '0');
setenv('i_rho', '0');
setenv('shaper', 'DPS_VSCI');

i_eps = str2double(getenv('i_eps'));
i_rho = str2double(getenv('i_rho'));
shaper = getenv('shaper');

ns = NetSetting('shaper',shaper);
[eps, rho] = ns.get_eps_rho(i_eps+1, i_rho+1);

pool_started = false;
while ~pool_started
    try
        delete(gcp('nocreate'));
        if strcmp(ns.where, 'cluster')
            pc = parcluster('local')
            parpool(pc, str2double(getenv('SLURM_CPUS_ON_NODE')))
        else
            parpool(2)
        end
        pool_started = true;
    catch
        continue;
    end
end

opt_res = ns.optim_shaper(eps, rho);

ns.i_eps = i_eps; ns.i_rho = i_rho;
res_prefix = strcat(shaper, '_opt');
fn = ns.get_fn(res_prefix, 'ind');
save(fn, '-struct', 'opt_res');