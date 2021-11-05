[mp] = rollingMP(data, w);

 h = heatmap(mp,'GridVisible','off');

Y = squareform(mp);
Z = linkage(Y);

dendrogram(Z);

T = cluster(Z,'cutoff',.999)';
T = cluster(Z,'maxclust',2)';


h = plot(1:length(data)-w+1, data(1:end-w+1));
cd = colormap('cool'); % take your pick (doc colormap)
cd = interp1(linspace(min(T),max(T),length(cd)),cd,T); % map color to y values
cd = uint8(cd'*255); % need a 4xN uint8 array
cd(4,:) = 255; % last column is transparency
set(h.Edge,'ColorBinding','interpolated','ColorData',cd);

hold on;
plot(length(data)-w+1:length(data), data(end-w+1:end),'color',[0.5,0.5,0.5]);
hold off;

drawnow;