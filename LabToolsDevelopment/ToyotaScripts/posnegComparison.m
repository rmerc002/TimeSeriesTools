figure;
tiledlayout(2,1);

ax1 = nexttile();
hold on;
plot(dataAnom(:,5));
plot(dataAnom(:,10));
plot(dataAnom(:,15));
hold off;
title('Anom');

ax2 = nexttile();
hold on;
plot(dataBase(:,5));
plot(dataBase(:,10));
plot(dataBase(:,15));
hold off;
title('Base');

linkaxes([ax1, ax2],'xy');



figure;
tiledlayout(2,1);

ax1 = nexttile();
hold on;
plot(posAnom);
plot(negAnom);
hold off;
title('Anom');
set(gca,'TickDir','out');
box off;

ax2 = nexttile();
hold on;
plot(posBase);
plot(negBase);
hold off;
title('Base');
set(gca,'TickDir','out');
box off;

linkaxes([ax1, ax2],'xy');