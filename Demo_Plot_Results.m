clear all
clc
close all
startup; % import EEGLab 
%% plot the figures: original ones
load results.mat
out = outt{1};
figure('visible','on')
set(gcf,'outerposition',get(0,'screensize'))
load chanlocs64.mat;
timeset = 500;
% HC
for mod = 1:4
    for comp = 1:3
        subplot(4,7,(mod-1)*7+comp);
        if mod == 1
            topoplot(U0{1}{mod}(:,comp),chanlocs64,'maplimits',[0  max(U0{1}{mod}(:,comp))]);
            %colorbar;
            %caxis([min(U0{1}{mod}(:,comp)),max(U0{1}{mod}(:,comp))]) 
        elseif mod == 2
            fl=4;fh=30;
            fIndex=linspace(fl,fh,130);
            plot(fIndex,U0{1}{mod}(:,comp),'color','b');
            set(gca,'XLim',[fl fh]);
            set(gca,'XTick',[4 10 15 20 25 30]);
            set(gca,'XTickLabel',[4 10 15 20 25 30])
%             xlabel('Frquency/Hz','fontsize',14);
        elseif mod == 3
            plot(1:timeset, U0{1}{mod}(:,comp),'b'); 
            xlim([1 timeset]);
        elseif mod == 4
            scatter(1:19,U0{1}{mod}(:,comp),'b');
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[1 19]);
            set(gca,'XTick',[1 5 10 15 19]);
            set(gca,'XTickLabel',[1 5 10 15 19])
        end
    end
end
% MDD
for mod = 1:4
    for comp = 1:4
        subplot(4,7,(mod-1)*7+comp+3);
        if mod == 1
            topoplot(U0{2}{mod}(:,comp),chanlocs64,'maplimits',[0  max(U0{2}{mod}(:,comp))]);
            %colorbar;
            %caxis([min(U0{2}{mod}(:,comp)),max(U0{2}{mod}(:,comp))]) 
        elseif mod == 2
            fIndex=linspace(4,30,130);
            plot(fIndex,U0{2}{mod}(:,comp),'color','b');
            set(gca,'XLim',[fl fh]);
            set(gca,'XTick',[4 10 15 20 25 30]);
            set(gca,'XTickLabel',[4 10 15 20 25 30])
%             xlabel('Frquency/Hz','fontsize',14);
        elseif mod == 3
            plot(1:timeset, U0{2}{mod}(:,comp),'b'); 
            xlim([1 timeset]);
        elseif mod == 4
            scatter(1:20,U0{2}{mod}(:,comp),'b');
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[1 20]);
            set(gca,'XTick',[1 5 10 15 20]);
            set(gca,'XTickLabel',[1 5 10 15 20])
        end
    end
end
%% plot the figures: recovered ones
[FMS1, ~, ~, perm] = score(ktensor(ones(3,1),U0{1}),ktensor(ones(3,1),out.U{1}),'lambda_penalty',false);
figure('visible','on')
set(gcf,'outerposition',get(0,'screensize'))
load chanlocs64.mat;
%HC
for mod = 1:4
    ps = 0;
    for comp = perm
        ps = ps + 1;
        subplot(4,7,(mod-1)*7+ps);
        %subplot(8,7,(mod-1+4)*7+comp);
        if mod == 1       
            topoplot(out.U{1}{mod}(:,comp),chanlocs64,'maplimits',[0  max(out.U{1}{mod}(:,comp))]);
            %colorbar;
            %caxis([min(out.U{1}{mod}(:,comp)),max(out.U{1}{mod}(:,comp))]) 
        elseif mod == 2
            fIndex=linspace(4,30,130);
            plot(fIndex,out.U{1}{mod}(:,comp),'color','b');
            set(gca,'XLim',[fl fh]);
            set(gca,'XTick',[4 10 15 20 25 30]);
            set(gca,'XTickLabel',[4 10 15 20 25 30])
%             xlabel('Frquency/Hz','fontsize',14);
        elseif mod == 3
            plot(1:timeset, out.U{1}{mod}(:,comp),'b'); 
            xlim([1 timeset]);
        elseif mod == 4
            scatter(1:19,out.U{1}{mod}(:,comp),'b');
            ymax = max(out.U{1}{mod}(:,comp));
            if ymax < 1,ymax = 1;end
            set(gca,'YLim',[0 ymax]);
            set(gca,'XLim',[1 20]);
            set(gca,'XTick',[1 5 10 15 19]);
            set(gca,'XTickLabel',[1 5 10 15 19])
        end
    end
end
% MDD
[FMS2, ~, ~, perm] = score(ktensor(ones(4,1),U0{2}),ktensor(ones(4,1),out.U{2}),'lambda_penalty',false);% post-permuation
for mod = 1:4
    ps = 0;
    for comp = perm
        ps = ps + 1;
        subplot(4,7,(mod-1)*7+ps+3);
        %subplot(8,7,(mod-1+4)*7+comp+3);
        if mod == 1
            topoplot(out.U{2}{mod}(:,comp),chanlocs64,'maplimits',[0  max(out.U{2}{mod}(:,comp))]);
            %colorbar;
            %caxis([min(out.U{2}{mod}(:,comp)),max(out.U{2}{mod}(:,comp))]) 
        elseif mod == 2
            fIndex=linspace(4,30,130);
            plot(fIndex,out.U{2}{mod}(:,comp),'color','b');
            set(gca,'XLim',[fl fh]);
            set(gca,'XTick',[4 10 15 20 25 30]);
            set(gca,'XTickLabel',[4 10 15 20 25 30])
%             xlabel('Frquency/Hz','fontsize',14);
        elseif mod == 3
            plot(1:timeset, out.U{2}{mod}(:,comp),'b');
            xlim([1 timeset]);
        elseif mod == 4
            scatter(1:20,out.U{2}{mod}(:,comp),'b');
            ymax = max(out.U{2}{mod}(:,comp));
            if ymax < 1,ymax = 1;end
            set(gca,'YLim',[0 ymax]);
            set(gca,'XLim',[1 20]);
            set(gca,'XTick',[1 5 10 15 20]);
            set(gca,'XTickLabel',[1 5 10 15 20])
        end
    end
end