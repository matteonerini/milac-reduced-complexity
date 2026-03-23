clear; clc;
rng(3);
tic;

NRs = 32:32:32*8;
NTs = 32:32:32*8;
SNR_dB = 10; % 0 or 10
nMonte = 1e1;
NSs = [4,8,12,16];

SNR = db2pow(SNR_dB);
Y0 = 1 / 50;

C = nan(nMonte,length(NRs),length(NSs));
R_digi = nan(nMonte,length(NRs),length(NSs));
R_full = nan(nMonte,length(NRs),length(NSs));
R_stem = nan(nMonte,length(NRs),length(NSs));

for iMonte = 1:nMonte
    for iNR = 1:length(NRs)
        NR = NRs(iNR);
        NT = NTs(iNR);

        H = sqrt(1/2) * (randn(NR,NT) + 1i * randn(NR,NT));

        [U_tmp,Sigma,V_tmp] = svd(H);
        alpha = 2 * pi * rand;
        U = U_tmp * exp(1i * alpha);
        V = V_tmp * exp(1i * alpha);
        sigma = diag(Sigma)';
        lambda = sigma.^2;

        for iNS = 1:length(NSs)
            NS = NSs(iNS);

            % Capacity
            p_NS = waterfill(1,4./(SNR*lambda(1:NS)));
            C(iMonte,iNR,iNS) = sum(log2(1 + SNR .* p_NS .* lambda(1:NS) ./ 4));

            % Digital
            WWH_NS = V_tmp(:,1:NS) * diag(p_NS) * V_tmp(:,1:NS)';
            R_digi(iMonte,iNR,iNS) = real(log2(det(eye(NR) + SNR * (H * WWH_NS * H') / 4)));

            % Fully-connected MiLAC
            BF_full = func_opt_fully_tx(V,NS);
            F_inv_full = inv(1i * BF_full / Y0 + eye(NS+NT));
            F_full = F_inv_full(NS+1:NS+NT,1:NS);

            BG_full = func_opt_fully_rx(U,NS);
            G_inv_full = inv(1i * BG_full / Y0 + eye(NR+NS));
            G_full = G_inv_full(NR+1:NR+NS,1:NR);
            
            pow_full = p_NS' .* abs(diag(G_full*H*F_full)).^2;
            pow_and_int_full = sum(abs(G_full*H*F_full*diag(sqrt(p_NS))).^2,2);
            G_norm_full = vecnorm(G_full,2,2);
            SINR_full = pow_full ./ (pow_and_int_full - pow_full + G_norm_full .^2 / SNR);
            R_full(iMonte,iNR,iNS) = sum(log2(1 + SINR_full));

            % Stem-connected MiLAC
            BF_stem = func_opt_stem_tx(V(:,1:NS));
            F_inv_stem = inv(1i * BF_stem / Y0 + eye(NS+NT));
            F_stem = F_inv_stem(NS+1:NS+NT,1:NS);

            BG_stem = func_opt_stem_rx(U(:,1:NS));
            G_inv_stem = inv(1i * BG_stem / Y0 + eye(NR+NS));
            G_stem = G_inv_stem(NR+1:NR+NS,1:NR);
            
            pow_stem = p_NS' .* abs(diag(G_stem*H*F_stem)).^2;
            pow_and_int_stem = sum(abs(G_stem*H*F_stem*diag(sqrt(p_NS))).^2,2);
            G_norm_stem = vecnorm(G_stem,2,2);
            SINR_stem = pow_stem ./ (pow_and_int_stem - pow_stem + G_norm_stem .^2 / SNR);
            R_stem(iMonte,iNR,iNS) = sum(log2(1 + SINR_stem));

        end
     
    end
end

toc;

C_mean = squeeze(mean(C));
R_digi_mean = squeeze(mean(R_digi));
R_full_mean = squeeze(mean(R_full));
R_stem_mean = squeeze(mean(R_stem));

%% Plot
figure('DefaultAxesFontSize',12);
LineW = 1.8; MarkS = 8;
hold on;
plot(NRs, R_stem_mean(:,4),'-','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC');
plot(NRs, R_stem_mean(:,3),'-','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC');
plot(NRs, R_stem_mean(:,2),'-','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC');
plot(NRs, R_stem_mean(:,1),'-','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC');
plot(NRs, R_full_mean(:,4),'--','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC');
plot(NRs, R_full_mean(:,3),'--','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC');
plot(NRs, R_full_mean(:,2),'--','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC');
plot(NRs, R_full_mean(:,1),'--','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC');
plot(NRs, C_mean(:,4),'^k','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','{\itC}, 16 streams');
plot(NRs, C_mean(:,3),'ok','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','{\itC}, 12 streams');
plot(NRs, C_mean(:,2),'sk','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','{\itC}, 8 streams');
plot(NRs, C_mean(:,1),'xk','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','{\itC}, 4 streams');
grid on; box on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('Number of antennas');
ylabel('Achievable rate [bps/Hz]');
legend('Location','northwest','NumColumns',3,'FontSize',10);
xlim([min(NRs) max(NRs)])
xticks(NRs)
if SNR_dB == 0
    ylim([10 70])
elseif SNR_dB == 10
    ylim([20 120])
end
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);