clear
clc
close all


%% a path-loss simulation of 3GPP UMa NLOS

% Constant
umi_los_sf_std_3gpp =6; % 3GPP TR 38.901 shadow fading standard deviation.
h_bs = 25; % Transmitter height 
h_ue = 1.5; % User device height
c = 3e8; % light speed
monte_num_trail = 1000;

% path-loss
d2D = 10:10:5000;
fc = [0.5e9, 6e9, 28e9, 47e9, 60e9];
pl_uma_nlos_3gpp = zeros(length(fc), length(d2D));
for fc_idx = 1:length(fc)
    curr_fc = fc(fc_idx); 
    for d2D_idx = 1:length(d2D)
        curr_d2D = d2D(d2D_idx);  
        curr_d3D = sqrt(curr_d2D^2 + (h_bs-h_ue)^2);
        d_bp = 4*(h_bs-1)*(h_ue-1)*curr_fc/c;
        for monte_idx = 1:monte_num_trail 
            pl_uma_nlos_3gpp_monte(fc_idx, d2D_idx) = 35.3 + 22.4 + 21.3*log10(curr_fc) + 20*log10(curr_d3D) - 0.3*(h_ue - 1.5);
        end
        pl_uma_nlos_3gpp(fc_idx,d2D_idx) = mean(pl_uma_nlos_3gpp_monte(fc_idx,d2D_idx,:)); 
    end
end

%% figure
for fc_idx = 1:length(fc)
    semilogx(d2D, pl_uma_nlos_3gpp(fc_idx,:)); hold on;
end
legend('3GPP UMa NLOS, 0.5GHz', '3GPP UMa NLOS, 6GHz', '3GPP UMa NLOS, 28GHz', '3GPP UMa NLOS, 47GHz', '3GPP UMa NLOS, 60GHz');
xlabel('T-R separation (m)'); ylabel('Path-loss (dB)')
