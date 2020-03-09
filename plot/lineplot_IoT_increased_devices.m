clc; clear; close;

subfolder = 'IoT/postproc';
res_prefix_dps_vsci = 'DPS_VSCI_iid';
res_prefix_vsci_csci = 'VSCI_CSCI_iid';

figure;
[wd, ht] = deal(4, 3);
set(gcf, 'PaperPosition', [0 0 wd ht]); %Position plot at left hand corner with width 6 and length 18.
set(gcf, 'PaperSize', [wd ht]);
ax_range = [0 1 1e-4 1e4];

%% IID input
ns = NetSetting('synthetic',false,'device','Sleep','where','local',...
                'var_out_interval',false);

res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));
% ns.coeff_var = ns.get_cv(ns.in_sizes(2:end), ns.in_rates(2:end));

qsz_actual_all_dps_vsci = res_dps_vsci.delay_all_DPS;
qsz_actual_all_vsci = res_vsci_csci.delay_all_VSCI;
qsz_actual_all_csci = res_vsci_csci.delay_all_CSCI;

subplot(1,3,1); 
myPlot('lineplot',ns,qsz_actual_all_vsci,' VSCI',qsz_actual_all_csci,' CSCI*',qsz_actual_all_dps_vsci, 'delay', ax_range);

%%
ns = NetSetting('synthetic',false,'device','Cam','where','local',...
                'var_out_interval',false);

res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));
% ns.coeff_var = ns.get_cv(ns.in_sizes(2:end), ns.in_rates(2:end));

qsz_actual_all_dps_vsci = res_dps_vsci.delay_all_DPS;
qsz_actual_all_vsci = res_vsci_csci.delay_all_VSCI;
qsz_actual_all_csci = res_vsci_csci.delay_all_CSCI;

subplot(1,3,2); 
myPlot('lineplot',ns,qsz_actual_all_vsci,' VSCI',qsz_actual_all_csci,' CSCI*',qsz_actual_all_dps_vsci, 'delay', ax_range);

%%
ns = NetSetting('synthetic',false,'device','SleepCamMerged','where','local',...
                'var_out_interval',false);
res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));
% ns.coeff_var = ns.get_cv(ns.in_sizes(2:end), ns.in_rates(2:end));

qsz_actual_all_dps_vsci = res_dps_vsci.delay_all_DPS;
qsz_actual_all_vsci = res_vsci_csci.delay_all_VSCI;
qsz_actual_all_csci = res_vsci_csci.delay_all_CSCI;

subplot(1,3,3); 
myPlot('lineplot',ns,qsz_actual_all_vsci,'VSCI',qsz_actual_all_csci,'CSCI*',qsz_actual_all_dps_vsci, 'delay', ax_range);

