vid1 = VideoReader("C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\MouseGait\Data\Coordinate Data-20221214T184034Z-001\Coordinate Data\Diseased\FILE1086_D194_C3_aMSH-20rpm-5min_Trim.mp4");
vid2 = VideoReader("C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\MouseGait\Data\Coordinate Data-20221214T184034Z-001\Coordinate Data\Diseased\FILE1087_D194_A1_aMSH-20rpm-5min_Trim.mp4");

ts = zeros(30,1);
ts(11:20) = ones(10,1);
ts = ts + randn(30,1)*1e-1;
subsequenceIndex = 11;
subsequenceLength = 10;
contextLength = 10;
sswc1 = subsequence(ts, subsequenceIndex, subsequenceLength, contextLength);
sswc2 = subsequence(ts, subsequenceIndex, subsequenceLength, contextLength);

numFrames = sswc1.length;%vid1.NumberOfFrames
startIndexVideo = 200;
endIndexVideo = startIndexVideo + numFrames-1;
videofig(numFrames, @(frm) redraw(frm, startIndexVideo, endIndexVideo, vid1, vid2, sswc1, sswc2));
redraw(1, startIndexVideo, endIndexVideo, vid1, vid2, sswc1, sswc2);

% v1 = VideoReader("C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\MouseGait\Data\Coordinate Data-20221214T184034Z-001\Coordinate Data\Diseased\FILE1086_D194_C3_aMSH-20rpm-5min_Trim.mp4")
% v2 = VideoReader("C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\MouseGait\Data\Coordinate Data-20221214T184034Z-001\Coordinate Data\Diseased\FILE1087_D194_A1_aMSH-20rpm-5min_Trim.mp4")
% i1 = 0;
% i2 = 0;
% fig = gcf;
% ax1 = subplot(2,2,1, 'Parent', fig);
% ax2 = subplot(2,2,2, 'Parent', fig);
% while i1 < v1.NumberOfFrames && i2 < v2.NumberOfFrames
%         if i1 < v1.NumberOfFrames
%             i1 = i1+1;
%             if ishandle(ax1)
%               image(ax1, v1.read(i1));
%             else
%               break;    %axes is gone, figure is probably gone too
%            end
%         end
%         if i2 < v2.NumberOfFrames
%             i2 = i2+1;
%             if ishandle(ax2)
%               image(ax2, v2.read(i2));
%             else
%               break;    %axes is gone, figure is probably gone too
%             end
%         end
%     drawnow
%     end