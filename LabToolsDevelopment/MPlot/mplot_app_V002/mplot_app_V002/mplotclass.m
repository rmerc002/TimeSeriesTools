classdef mplotclass
   properties
      m = 100;
      cc = 0.9;
      fig;
      ax_mplot;
      ax_A;
      ax_B;
      ax_tsA;
      ax_tsB;
      data;
      ds_data;
      mplot;
      mslider;
      ccslider;
      sw;
      paa_factor = 1;
      dcm;
      dcm_enabled = false;
      
   end
   methods
       function startup(obj, timeseries)
            clc;
            obj.data = timeseries;
            obj.ds_data = obj.data;
            obj.fig = uifigure('Position',[100 100 1000 1000]);
            pos2 = [35, 10, 700, 700];
            obj.ax_mplot = uiaxes(obj.fig,'Position', pos2);
            obj.ax_mplot.FontSize = 11;
            obj.ax_tsA = uiaxes(obj.fig,'Position',[750,450,220,80]);
            obj.ax_tsB = uiaxes(obj.fig,'Position',[750,350,220,80]);
            obj.ax_tsA.YTick = [];
            obj.ax_tsB.YTick = [];
            set(obj.ax_tsA,'visible','off');
            set(obj.ax_tsB,'visible','off');
            
            obj.mslider = uislider(obj.fig,...
            'Position',[750,700,200,23],...
            'MajorTicks',[5,50,100,150,200,250,300],...
            'ValueChangedFcn',@(src,event) sliderCallback);
            obj.mslider.Limits = [5 300];
            obj.mslider.Value = 100;
        
             obj.ccslider = uislider(obj.fig,...
            'Position',[750,650, 200,23],...
            'MajorTicks',[0,0.2,0.4,0.6,0.8,1],...
            'ValueChangedFcn',@(src,event) ccsliderCallback);
            obj.ccslider.Limits = [0 1];
            obj.ccslider.Value = 0.9;
                        
            uibutton(obj.fig,'push',...
               'Position',[750,570,200,20],...
               'Text', 'Explore Subsequences (ON/OFF)', ...
               'ButtonPushedFcn', @(btn,event) subButtonPushed);
           
            obj.dcm = datacursormode(obj.fig); 
            
           
            function sliderCallback(src,event)
                obj.m = round(obj.mslider.Value);
                plotMplot(obj)
            end

            function ccsliderCallback(src,event)
                obj.cc = obj.ccslider.Value;
                plotMplot(obj)
            end
            
            plotMplot(obj);
           
           function subButtonPushed
               if obj.dcm_enabled
                   obj.dcm_enabled = false;
                   obj.dcm.Enable = 'off';
                   cla(obj.ax_tsA);
                   cla(obj.ax_tsB);
                   set(obj.ax_tsA,'visible','off');
                   set(obj.ax_tsB,'visible','off');
               else
                   obj.dcm_enabled = true;
                   obj.dcm.Enable = 'on';
                   obj.dcm.UpdateFcn = @displayCoordinates;
               end
           end
           
           function txt = displayCoordinates(~,info)
                x = info.Position(1);
                y = info.Position(2);
                txt = ['(' num2str(x) ', ' num2str(y) ')'];
                exploreSubsequences(obj, x, y);
           end
            
       end
       function plotMplot(obj)
            [obj.mplot, obj.paa_factor, ~] = SPLAT(obj.data, obj.m, nan, 0, 1,0);
            if obj.paa_factor > 1
                data_newlength = floor(length(obj.data)/obj.paa_factor);
                obj.ds_data = paa(obj.data, data_newlength);
            else
                obj.ds_data = obj.data;
            end
            
            pos1 = [35+25, 12 + 700, 700 - 25, 30];
            
            pos3 = [10, 10 + 5, 30, 700 - 5];

            
            histmaximum = max(obj.mplot,[],'all');
            obj.mplot = obj.mplot > obj.cc *histmaximum;
            
            imagesc(obj.mplot, 'Parent', obj.ax_mplot)
            colormap(obj.ax_mplot,gray)
            xt = get(obj.ax_mplot, 'XTick');
            set(obj.ax_mplot, 'XTick',[min(xt):floor((max(xt)-min(xt))/10):max(xt)]);
            yt = get(obj.ax_mplot, 'YTick');
            set(obj.ax_mplot, 'YTick',[min(yt):floor((max(yt)-min(yt))/10):max(yt)]);
            %obj.ax_mplot.DataAspectRatio = [1,1,1];
            
            
            obj.ax_A = uiaxes(obj.fig,'Position', pos1);
            plot(obj.ds_data, 'Parent', obj.ax_A);
            obj.ax_A.XTick = [];
            obj.ax_A.YTick = [];
            obj.ax_A.Color = 'None';
            set(obj.ax_A, 'visible', 'off')
            
            
            obj.ax_B = uiaxes(obj.fig,'Position', pos3);
            plot(obj.ds_data, 'Parent', obj.ax_B);
            obj.ax_B.XTick = [];
            obj.ax_B.YTick = [];
            obj.ax_B.Color = 'None';
            view(obj.ax_B,[90,-90])
            set(obj.ax_B,'xdir','reverse','ydir','reverse', 'XAxisLocation','top');
            set(obj.ax_B, 'visible', 'off')
            %linkaxes([obj.ax_mplot, obj.ax_A],'x');
       end
       
       function exploreSubsequences(obj, x, y)
            p = obj.paa_factor;
            subA = obj.data(max(1,x*p-obj.m):min(x*p+obj.m,end));
            subB = obj.data(max(1,y*p-obj.m):min(y*p+obj.m,end));
            set(obj.ax_tsA,'visible','on', 'Color', 'none');
            set(obj.ax_tsB,'visible','on', 'Color', 'none');
            plot(obj.ax_tsA, subA)
            title(obj.ax_tsA, sprintf('T_A^{%.0f}',x*p));
            plot(obj.ax_tsB, subB)
            title(obj.ax_tsB, sprintf('T_B^{%.0f}',y*p));
            xt = get(obj.ax_tsA, 'XTick');
            set(obj.ax_tsA, 'XTick',[min(xt) max(xt)], 'XTickLabel',{'x-m/2','x+m/2'});
            yt = get(obj.ax_tsB, 'XTick');
            set(obj.ax_tsB, 'XTick',[min(yt) max(yt)], 'XTickLabel',{'y-m/2','y+m/2'});
       end
   end
end