function [x, y, z] = rfPlot(params, RF, parent, peak)
% rfPlot - script to visualize cropped RF
%
% [x, y, z] = rfPlot(params, RF, [parent=axes in new figure]);
%
% Will produce a plot illustrating the estimated location of a 
% 2D Gaussian receptive field, with a grid.
%
% RF is a 
% the 'parent' argument directs where to display the plot. Default
% is to create a new figure with its own axes.
%
% 2006/02 SOD: wrote it.
% 2006/09 RAS: added optional 'parent' argument, so you can
% direct the plot to a subplot axes.
% 2008/06 RAS: updated calculation of RF grid to use the X, Y sample points
% in params.analysis.X (and .Y). This replaces a previous method using the
% sample rate; I've found that when you recompute the stimulus (e.g.
% rmRecomputeParams), X/Y sample points for which no stimulus was presented
% are omitted from the analysis. We use this sampling grid to be
% consistent, and prevent bugs in code like rmPlotGUI.
% 2008/07 SOD: reverted back to original. See comments below at the 
% relevant code. Must validate rmRecomputeParams.
if ~exist('parent','var') || isempty(parent), figure; parent = gca;      end;
if ~exist('peak','var'), peak = [];      end;

if params.analysis.numberStimulusGridPoints==441 || params.analysis.numberStimulusGridPoints==987 || params.analysis.numberStimulusGridPoints==49
    [x,y]=meshgrid(0:0.01:max(params.analysis.X), 0:0.01:max(params.analysis.Y));
    z    = zeros(size(x));
    %z(params.stim(1).instimwindow) 
    z = RF;
    z    = reshape(z,size(x));
%     if params.analysis.X(1)==0.05
%         [newx, newy]=meshgrid(min(params.analysis.X):0.05:(0.05+max(params.analysis.X)), min(params.analysis.Y):0.05:(0.05+max(params.analysis.Y)));
%         newz    = zeros(size(newx));
%         [tmp,a]=intersectCols(round([newx(:) newy(:)]'.*20), round([x(:) y(:)]'.*20));
%         newz(a)=z;
%         x=newx;
%         y=newy;
%         z=newz;
%     else


%         newx=zeros(size(x,1)+1, size(x,2)+1);
%         newy=newx;
%         newz=newx;
%         newx(1:end-1, 1:end-1)=x;
%         newy(1:end-1, 1:end-1)=y;
%         newx(:, end)=newx(:,end-1)+0.1;
%         newx(end,:)=newx(end-1,:);
%         newy(end,:)=newy(end-1,:)+0.1;
%         newy(:,end)=newy(:,end-1);
%         newy(1:end-1, 1:end-1)=y;        
%         newz(1:end-1, 1:end-1)=z;
%         x=newx;
%         y=newy;
%         z=newz;


%    end
else
    [x,y] = prfSamplingGrid(params);
    z    = nan(size(x));
    z(params.stim(1).instimwindow) = RF;
    z    = reshape(z,size(x));
end

% we need to flip the RF vertically because positive y-values correspond
% to upper visual field (and hence should be plotted high), but high vales
% in an image matrix are plotted low
%z = flipud(z);

%z=z.*y;
z=z-max(z(:));

% plot
axes(parent);
cla;
hold on;
surf(x,y,z,'LineStyle','none');% colormap(gray)

%draw lines at every degree
mylines = 0; %[-floor(params.analysis.fieldSize):floor(params.analysis.fieldSize)];
for ll = 1:numel(mylines),
    for n=1:2,
        if n==1,
            ii = find(x==mylines(ll) & isfinite(z));
            [xs, is] = sort(y(ii));
        else
            ii = find(y==mylines(ll) & isfinite(z));
            [xs, is] = sort(x(ii));
        end;
        if ~isempty(is),
            ii = ii(is);
            h  = line(x(ii),y(ii),z(ii));
            if mylines(ll)==0,
                set(h,'LineWidth',1,'Color',[0 0 0]);
            else
                set(h,'LineWidth',0.5,'Color',[0 0 0]);
            end;
        end;
    end;
end;

minz = min(z(:));
maxz = peak;
if isempty(maxz), maxz = max(z(:)); end;
if isnan(maxz), maxz = 0.1; minz = -0.1; end;
if minz==maxz,
    minz = minz - 0.1;
    maxz = maxz + 0.1;
end;


% scale axis
axis([min(x(:)) max(x(:)) min(y(:)) max(y(:)) minz maxz]);
axis image;
colormap gray;

% scale colorbar to be centered on zero
caxis([0 1].*max(abs(minz),abs(maxz)));
caxis([min(z(:)) max(z(:))]);
hold on; plot(params.analysis.X, params.analysis.Y, 'r.')
[ymax, xmax]=ind2sub(size(z), find(z==max(z(:))));
hold on; plot(x(1,xmax), y(ymax,1), 'b.');

% axis labels
xlabel('x (deg)'); 
ylabel('y (deg)'); 
zlabel('BOLD amplitude (%/deg^2/sec)'); 
title('pRF profile');
hold off;

return
