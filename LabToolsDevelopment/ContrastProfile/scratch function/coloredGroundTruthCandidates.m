x = zeros(10,1000);
y = zeros(size(x));

figure;
plot(0,0);
hold on;
for i = 1:size(x,1)
    x(i,:) = 1:size(x,2);
    y(i,:) = sin(linspace(0,10*2*pi,size(x,2)));
    
    tempMin = min(y(i,:));
    tempMax = max(y(i,:));
    tempRange = max(1e-5, (tempMax-tempMin));

    cmapping = false(1,length(x));
    gti = randi(size(x,2));
    fprintf('ground truth index: %d\n',gti);
    cmapping(gti:end) = true;

    xfalse = x(i,:);
    xtrue = x(i,:);

    xfalse(cmapping) = nan;
    xtrue(~cmapping) = nan;

    yfalse = y(i,:);
    ytrue = y(i,:);

    yfalse(cmapping) = nan;
    ytrue(~cmapping) = nan;

    

    % scatter(x(cmapping==false),y(cmapping==false),3,'MarkerFaceColor',[0.75,0.75,0.75]);
    % scatter(x(cmapping==true),y(cmapping==true),3,'MarkerFaceColor',[0,0.2,0.8]);
    plot(xfalse,1-i+0.95*(yfalse-tempMin)/tempRange,'Color',[0.5,0.5,0.5]);
    plot(xtrue,1-i+0.95*(ytrue-tempMin)/tempRange,'Color',[0,0.8,0.2]);
    scatter(gti,1-i+0.95*(ytrue(gti)-tempMin)/tempRange,10,'MarkerFaceColor',[0,0.8,0.2], 'MarkerEdgeAlpha',0);
    
    
end
r = rectangle();
r.Position = [1, 1-i, ceil(size(x,2)/4), 5];
r.FaceColor = [1,1,1, 0.45];
r.LineStyle = '--';
r.EdgeColor = [0,0,0, 0.45];

r = rectangle();
r.Position = [ceil(size(x,2)/4)*3, 1-i, ceil(size(x,2)/4), 5];
r.FaceColor = [1,1,1, 0.45];
r.LineStyle = '--';
r.EdgeColor = [0,0,0, 0.45];
hold off;
