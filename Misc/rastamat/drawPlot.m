function  drawPlot(MFCC_Coef,d, hopTime,duration)

coefSR = 1/hopTime;
plotNum = 15;
h = zeros(plotNum,1);

figure('Name','MFFC Coefficients','NumberTitle','off');
s = strcat('Duration = ',num2str(duration),' s');
mTextBox = uicontrol('style','text');
set(mTextBox,'String',s,'Position',[0 400 300 20],'Units','normalized','FontUnits','normalized','FontSize',0.7);
h(1) = subplot(plotNum,1,1);
plot(d);

[k l] = size(d);
[m n] = size(MFCC_Coef(1,:));
d = downsample(d,ceil(k/(n)));

[m , n] = size(MFCC_Coef());
nCoeff = zeros(m,n);
ymax = -inf;
ymin = inf;
for i=1:1:m
    nCoeff(i,:) = zscore(MFCC_Coef(i,:));
    maxValu = max(nCoeff(i,:)); minValue = min(nCoeff(i,:));
    if(maxValu > ymax)
       ymax = maxValu;
    end
    if(minValue < ymin)
        ymin = minValue;
    end
end

h(2)=subplot(plotNum,1,3);
plot(d);
set(gca, 'box', 'off');

mTextBox2 = uicontrol('style','text');
set(mTextBox2,'String','Downsample Sound','Position',[1 324 50 20],'Units','normalized','FontUnits','normalized','FontSize',0.45);
mTextBox3 = uicontrol('style','text');
set(mTextBox3,'String','Original Sound','Position',[1 365 50 20],'Units','normalized','FontUnits','normalized','FontSize',0.45);
mTextBox1 = uicontrol('style','text');
s = strcat('SR of Coefs = ',num2str(coefSR),' HZ');
set(mTextBox1,'String',s,'Position',[200 400 300 20],'Units','normalized','FontUnits','normalized','FontSize',0.7);
set(gca,'XTick',[]);

startPlot = 3; % three plot has been draw
for(i=1:1:(plotNum-startPlot))
    h(i+startPlot) = subplot( plotNum,1,(i+startPlot) );
    plot(nCoeff(i,:));
    mTextBox = uicontrol('style','text');
    set(mTextBox,'String',num2str(i),'Position',[25 323-23.5*i 20 20],'Units','normalized','FontUnits','normalized','FontSize',0.60);
    ylim([ymin,ymax]);
    set(gca, 'box', 'off');
    if((i+startPlot)<plotNum)
        set(gca,'XTick',[]);
    end
end
linkaxes( h(2:plotNum), 'x' );
end
