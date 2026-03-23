clear; clc;

NTs = 32:32:32*8;
NSs = [4,8,12,16];

MiLAC_full = nan(length(NTs),length(NSs));
MiLAC_stem = nan(length(NTs),length(NSs));
for iNS = 1:length(NSs)
    NS = NSs(iNS);
    MiLAC_stem(:,iNS) = NS*(2*NTs+1);
    MiLAC_full(:,iNS) = (NS+NTs).*(NS+NTs+1)/2;
end

figure('defaultaxesfontsize',12)
LineW = 1.8; MarkS = 8;
plot(NTs,MiLAC_stem(:,4),'-^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC')
hold on;
plot(NTs,MiLAC_stem(:,3),'-o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC')
plot(NTs,MiLAC_stem(:,2),'-s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC')
plot(NTs,MiLAC_stem(:,1),'-x','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Stem-c. MiLAC')
set(gca,'ColorOrderIndex', 1)
plot(NTs,MiLAC_full(:,4),'--^','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC, 16 streams')
plot(NTs,MiLAC_full(:,3),'--o','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC, 12 streams')
plot(NTs,MiLAC_full(:,2),'--s','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC, 8 streams')
plot(NTs,MiLAC_full(:,1),'--x','linewidth',LineW,'MarkerSize',MarkS,'DisplayName','Fully-c. MiLAC, 4 streams')
grid on; box on;
set(gca,'GridLineStyle',':','GridAlpha',0.5,'LineWidth',1.2);
xlabel('Number of antennas');
ylabel('Circuit complexity');
legend('location','northwest','NumColumns',2,'FontSize',10);
ax = gca;
ax.XTick = 0:32:32*8;
ax.XLim = [32 32*8];
ax.YTick = 0:2000:16000;
ax.YLim = [0 16000];
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1.5);