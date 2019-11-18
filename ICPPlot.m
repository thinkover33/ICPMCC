%  ICPPlot(X, Y, C); plots 2 data sets. Works only for 2D and 3D data sets.
%
%   Input
%   ------------------
%   X           Reference point set matrix NxD;
%   Y           Current postions of GMM centroids;

function ICPPlot(X, Y, gcf)

if nargin<2, error('icp_plot.m error! Not enough input parameters.'); end;
if nargin<3
    gcf = figure;
end
[m, d]=size(Y);
if d>3, error('icp_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;
if d<2, error('icp_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;

nXY = max(X);
nX1 = ceil(nXY(1));
nXY = max(Y);
nX2 = ceil(nXY(1));
nXYmax = max(nX1, nX2);
%
figure(gcf);
% for 2D case
if d==2,
    plot( X(:,2), nXYmax-X(:,1),'b.',Y(:,2), nXYmax-Y(:,1),'g.');
    axis equal;
    %  axis off; axis([-1.5 2 -1.5 2]);
    %    plot(X(:,1), X(:,2),'r.', Y(:,1), Y(:,2),'b.');
    %   axis off; axis([-1.5 2 -1.5 2]);
    %    data=nXYmax-X(:,1);
    %    legend('data','nXYmax-Y(:,1)',2);
    
else
    % for 3D case
    plot3(X(:,1),X(:,2),X(:,3),'b.','MarkerSize',1);
    hold on
    plot3(Y(:,1),Y(:,2),Y(:,3),'r.','MarkerSize',1);
    axis equal;
    % title('X data (red). Y GMM centroids (blue)');set(gca,'CameraPosition',[15 -50 8]);
%     ptX = pointCloud(X);
%     ptY = pointCloud(Y);
%     pcshow(ptX);
%     hold on;
%     pcshow(ptY);
%     axis equal;
end

drawnow;