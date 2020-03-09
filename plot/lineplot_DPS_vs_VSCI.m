clc; clear; close;

subfolder = 'DPS_vs_VSCI/postproc/old';
res_prefix_dps_vsci = 'DPS_VSCI_iid';
res_prefix_vsci_csci = 'VSCI_CSCI_iid';

scale_all = [0.01, 1, 5];

figure;
[wd, ht] = deal(4, 3);
set(gcf, 'PaperPosition', [0 0 wd ht]); %Position plot at left hand corner with width 6 and length 18.
set(gcf, 'PaperSize', [wd ht]);
set(gca, 'LooseInset', get(gca,'TightInset'));
ax_range = [1e-2 1 log(10^(-3.5)) log(1e3)];

n_subplots = length(scale_all);

%% IID input
i = 3;
ns = NetSetting('synthetic',true, 'scale',scale_all(i), 'where','local',...
                'var_out_interval',false,'scaled_alphabet',false, 'increased_alphabet',false);

res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));

qsz_actual_all_dps_vsci = res_dps_vsci.qsz_actual_all_DPS;
qsz_actual_all_vsci = res_vsci_csci.qsz_actual_all_VSCI;
qsz_actual_all_csci = res_vsci_csci.qsz_actual_all_CSCI;

% subplot(2,n_subplots,i); 
myPlot('lineplot',ns,qsz_actual_all_vsci,' VSCI',qsz_actual_all_csci,' CSCI*',qsz_actual_all_dps_vsci, 'qsz', ax_range);

delay_all_dps_vsci = res_dps_vsci.delay_all_DPS;
delay_all_vsci = res_vsci_csci.delay_all_VSCI;
delay_all_csci = res_vsci_csci.delay_all_CSCI;
% subplot(2,n_subplots,i+n_subplots); 
% myPlot('lineplot',ns,delay_all_vsci,' VSCI',delay_all_csci,' CSCI*',delay_all_dps_vsci, 'delay', ax_range);
