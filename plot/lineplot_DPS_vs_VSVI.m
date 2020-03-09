% clc; clear;

subfolder = 'DPS_vs_VSVI/postproc/old';
res_prefix_dps_vsvi = 'DPS_VSVI_iid';
res_prefix_vsvi_csvi = 'VSVI_CSVI_iid';

scale_all = [0.01, 1, 5];

figure;
[wd, ht] = deal(4, 3);
set(gcf, 'PaperPosition', [0 0 wd ht]); %Position plot at left hand corner with width 6 and length 18.
set(gcf, 'PaperSize', [wd ht]);
ax_range = [1e-2 1 10^(-3.5) 1e4];

n_subplots = length(scale_all);

%% IID input
i = 1;
ns = NetSetting('synthetic',true, 'scale',scale_all(i), 'where','local','var_out_interval',true);    

res_dps_vsvi = load(ns.get_rn(subfolder, res_prefix_dps_vsvi));
res_vsvi_csvi = load(ns.get_rn(subfolder, res_prefix_vsvi_csvi));

qsz_actual_all_dps_vsvi = res_dps_vsvi.delay_all_DPS;
qsz_actual_all_vsvi = res_vsvi_csvi.delay_all_VSVI;
qsz_actual_all_csvi = res_vsvi_csvi.delay_all_CSVI;

% subplot(1,n_subplots,i); hold on; grid on;
myPlot('lineplot',ns,qsz_actual_all_vsvi,'VSVI',qsz_actual_all_csvi,'CSVI*',qsz_actual_all_dps_vsvi, 'delay', ax_range);

