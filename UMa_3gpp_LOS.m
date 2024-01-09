clear
clc
close all

%% a path-loss simulation of 3GPP UMa LOS

% Constant
umi_los_sf_std_3gpp =4; % 3GPP TR 38.901 shadow fading standard deviation.
h_bs = 25; % Transmitter height 
h_ue = 1.5; % User device height
c = 3e8; % light speed
monte_num_trail = 1000;

% path loss
d2D = 10:10:5000;
fc = [0.5e9, 6e9, 28e9, 47e9, 60e9];
pl_uma_los_3gpp = zeros(length(fc), length(d2D));
for fc_idx = 1:length(fc)
    curr_fc = fc(fc_idx); 
    for d2D_idx = 1:length(d2D)
        curr_d2D = d2D(d2D_idx); 
        curr_d3D = sqrt(curr_d2D^2 + (h_bs-h_ue)^2); 
        d_bp = 4*h_bs*h_ue*curr_fc/c; 
        for monte_idx = 1:monte_num_trail 
            if curr_d2D <= d_bp
                pl_uma_los_3gpp_monte(fc_idx,d2D_idx,monte_idx) = 28.0 + 22*log10(curr_d3D) + 20*log10(curr_fc/1e9) + uma_los_sf_std_3gpp*randn(1); %PL1
            else
                pl_uma_los_3gpp_monte(fc_idx,d2D_idx,monte_idx) = 28.0 + 40*log10(curr_d3D) + 20*log10(curr_fc/1e9) - 9*log10(d_bp^2 + (h_bs-h_ue)^2) + uma_los_sf_std_3gpp*randn(1); % PL2
            end
            
        end
        pl_uma_los_3gpp(fc_idx,d2D_idx) = mean(pl_uma_los_3gpp_monte(fc_idx,d2D_idx,:)); 
    end
end

%% figure： d2D VS path-loss
for fc_idx = 1:length(fc)
    semilogx(d2D, pl_uma_los_3gpp(fc_idx,:)); hold on;
end
legend('3GPP UMa LOS, 0.5GHz', '3GPP UMa LOS, 6GHz', '3GPP UMa LOS, 28GHz', '3GPP UMa LOS, 47GHz', '3GPP UMa LOS, 60GHz');
xlabel('T-R separation (m)'); ylabel('Path-loss (dB)')
