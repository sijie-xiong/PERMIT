function main()
clc; clear; close all;

subfolder = 'IoT/postproc';

ns = NetSetting('synthetic',false,'device','Cam');
n = 200;
k = 4;
[i_eps,i_rho] = deal(5,3);
ax_range = [0 n 0 max(ns.data.sizes)+30];

%% Original iid trace
iid_trace = randsample(ns.data.sizes,n,true,ns.data.pmf);
% VSCI/CSCI*
res_prefix_dps_vsci = 'DPS_VSCI_iid';
res_prefix_vsci_csci = 'VSCI_CSCI_iid';
res_dps_vsci = load(ns.get_rn(subfolder, res_prefix_dps_vsci));
res_vsci_csci = load(ns.get_rn(subfolder, res_prefix_vsci_csci));

iid_trace_vsci = zeros(size(iid_trace));
for i = 1:n
    mu_opt = res_vsci_csci.mu_opt_all(i_rho,:);
%     pkt_sz = iid_trace(i);
%     idx = find(ns.data.sizes == pkt_sz);
    iid_trace_vsci(i) = randsample(ns.data.sizes, 1, true, mu_opt);
end
iid_trace_csci = ones(size(iid_trace))*res_vsci_csci.exp_pkt_sz_all(i_rho);
% DPS-VSCI
iid_trace_dps = zeros(size(iid_trace));
for i = 1:n
    Q_opt = res_dps_vsci.Q_opt_all(:,:,i_eps,i_rho);
    pkt_sz = iid_trace(i);
    idx = find(ns.data.sizes == pkt_sz);
    iid_trace_dps(i) = randsample(ns.data.sizes, 1, true, Q_opt(idx,:));
end
% DPS-VSVI

%% Original bursty trace
bursty_trace = ns.data.trace(k*n:(k+1)*n);
% VSCI/CSCI*
bursty_trace_vsci = zeros(size(bursty_trace));
for i = 1:n
    mu_opt = res_vsci_csci.mu_opt_all(i_rho,:);
%     pkt_sz = iid_trace(i);
%     idx = find(ns.data.sizes == pkt_sz);
    bursty_trace_vsci(i) = randsample(ns.data.sizes, 1, true, mu_opt);
end
bursty_trace_csci = ones(size(bursty_trace))*res_vsci_csci.exp_pkt_sz_all(i_rho);
% DPS-VSCI
bursty_trace_dps = zeros(size(bursty_trace));
for i = 1:n
    Q_opt = res_dps_vsci.Q_opt_all(:,:,i_eps,i_rho);
    pkt_sz = bursty_trace(i);
    idx = find(ns.data.sizes == pkt_sz);
    bursty_trace_dps(i) = randsample(ns.data.sizes, 1, true, Q_opt(idx,:));
end
% DPS-VSVI

my_trace_plot(iid_trace, ax_range, 'Orig. i.i.d. trace');
my_trace_plot(iid_trace_vsci, ax_range, 'i.i.d. trace after VSCI');
my_trace_plot(iid_trace_csci, ax_range, 'i.i.d. trace after CSCI*');
my_trace_plot(iid_trace_dps, ax_range, 'i.i.d. trace after \epsilon/\epsilon/\epsilon-DPS');
my_trace_plot(bursty_trace, ax_range, 'Orig. bursty trace');
my_trace_plot(bursty_trace_vsci, ax_range, 'bursty trace after VSCI');
my_trace_plot(bursty_trace_csci, ax_range, 'bursty trace after CSCI*');
my_trace_plot(bursty_trace_dps, ax_range, 'bursty trace after \epsilon/\epsilon/\epsilon-DPS');

end

function my_trace_plot(trace, ax_range, title_str)
figure;
[wd, ht] = deal(5, 3);
set(gcf, 'PaperPosition', [0 0 wd ht]); %Position plot at left hand corner with width 6 and length 18.
set(gcf, 'PaperSize', [wd ht]);
fontsize = 20;
lindwidth = 3;
bar(trace,'linewidth',lindwidth,'facecolor','#0072BD');
set(gca,'FontSize',fontsize+5);
title(title_str,'FontSize',fontsize+5,'FontWeight','bold');
xlabel('Time Slot','FontSize',fontsize+5,'FontWeight','bold'); 
ylabel('Packet Size (Bytes)','FontSize',fontsize+5,'FontWeight','bold'); 
axis(ax_range);
% subplot(2,4,1); bar(iid_trace); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('Orig. i.i.d. trace'); axis(ax_range);
% subplot(2,4,2); bar(iid_trace_vsci); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('i.i.d. trace after VSCI'); axis(ax_range);
% subplot(2,4,3); bar(iid_trace_csci); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('i.i.d. trace after CSCI*'); axis(ax_range);
% subplot(2,4,4); bar(iid_trace_dps); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('i.i.d. trace after \epsilon/\epsilon/\epsilon-DPS'); axis(ax_range);
% subplot(2,4,5); bar(bursty_trace); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('Orig. bursty trace'); axis(ax_range);
% subplot(2,4,6); bar(bursty_trace_vsci); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('bursty trace after VSCI'); axis(ax_range);
% subplot(2,4,7); bar(bursty_trace_csci); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('bursty trace after CSCI*'); axis(ax_range);
% subplot(2,4,8); bar(bursty_trace_dps); 
% xlabel('Time Slot'); ylabel('Packet Size (Bytes)'); title('bursty trace after \epsilon/\epsilon/\epsilon-DPS'); axis(ax_range);
end
