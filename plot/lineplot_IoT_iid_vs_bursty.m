clc; clear; close all;

subfolder = 'IoT/postproc';

figure;
[wd, ht] = deal(4, 3);
set(gcf, 'PaperPosition', [0 0 wd ht]); %Position plot at left hand corner with width 6 and length 18.
set(gcf, 'PaperSize', [wd ht]);
ax_range = [0.1 1 1e-4 1e4];

% n_subplots = length(scale_all);
%% IID input
ns = NetSetting('synthetic',false,'device','Cam','where','local');

res_prefix_dps_vsci = 'DPS_VSCI_iid';
res_prefix_vsci_csci = 'VSCI_CSCI_iid';
res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));

% qsz_actual_all_dps_vsci = res_dps_vsci.qsz_est_all;
% qsz_actual_all_vsci = res_vsci_csci.qsz_est_all_VSCI;
% qsz_actual_all_csci = res_vsci_csci.qsz_est_all_CSCI;
qsz_actual_all_dps_vsci = res_dps_vsci.delay_all_DPS;
qsz_actual_all_vsci = res_vsci_csci.delay_all_VSCI;
qsz_actual_all_csci = res_vsci_csci.delay_all_CSCI;

subplot(1,2,1); 
myPlot('lineplot',ns,qsz_actual_all_vsci,' VSCI',qsz_actual_all_csci,' CSCI*',qsz_actual_all_dps_vsci, 'delay', ax_range);

ns.input_typ = 'bursty';
res_prefix_dps_vsci = 'DPS_VSCI_bursty';
res_prefix_vsci_csci = 'VSCI_CSCI_bursty';
res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));

% qsz_actual_all_dps_vsci = res_dps_vsci.qsz_actual_all_DPS;
% qsz_actual_all_vsci = res_vsci_csci.qsz_actual_all_VSCI;
% qsz_actual_all_csci = res_vsci_csci.qsz_actual_all_CSCI;
qsz_actual_all_dps_vsci = res_dps_vsci.delay_all_DPS;
qsz_actual_all_vsci = res_vsci_csci.delay_all_VSCI;
qsz_actual_all_csci = res_vsci_csci.delay_all_CSCI;
% qsz_actual_all_dps_vsci(end,end) = 1e-5;

subplot(1,2,2); 
myPlot('lineplot',ns,qsz_actual_all_vsci,' VSCI',qsz_actual_all_csci,' CSCI*',qsz_actual_all_dps_vsci, 'delay', ax_range);

