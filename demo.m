%   3D point sets registration
clc;
clear;
close all;
%%%%%%A*MData->DData找对应点

% load('dragon');
% MData = dragon{3,1};
% DData = dragon{2,1};

% load('bunny');
% MData = bunny{1,1};
% DData = bunny{3,1};

load('bunny');
MData = bunny{2,1};
DData = bunny{3,1};
% 参数：80*0.97

% [NMap,D] = size(Map0);
% Map0 = Map0(floor(NMap/3):floor(NMap/3*2),:);
% MData = Pcd0;
% DData = Map0;
% clear Map0 Pcd0
%%%%%%%%仿真部分%%%%%%%%%%%
% [M,N] = size(DData);
% DData = DData(1:M*0.8,:);%删除一些点
%%%%%%%%%画出形状%%%%%%%
figure
plot3(DData(:,1),DData(:,2),DData(:,3),'r.')
axis equal;
grid on
xlabel('x');
ylabel('y');
zlabel('z');
figure
plot3(MData(:,1),MData(:,2),MData(:,3),'b.')
axis equal;
grid on
xlabel('x');
ylabel('y');
zlabel('z');

%%%%%%%%%%%%%%%%%%%%%%%% ICP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp('ICP');
% [R1, t1, TCorr1, MSE1, TData1] = ICP(MData, DData);
% figure
% plot3(DData(:,1),DData(:,2),DData(:,3),'r.')
% hold on
% plot3(TData1(:,1),TData1(:,2),TData1(:,3),'b.')
% axis equal;
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('ICP result');

[R2, t2, TCorr2, MSE2, TData2] = ICP_MCC(MData, DData);
figure
plot3(DData(:,1),DData(:,2),DData(:,3),'r.') 
hold on
plot3(TData2(:,1),TData2(:,2),TData2(:,3),'b.')
axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('ICP MCC result')

[R3, t3, TCorr3, MSE3, TData3] = ICP_MCC_pt2pl(MData, DData);
ICPPlot(TData3,DData); %YData是不动的一方，所以与原图一样
axis equal;
axis off
title('ICPMCC p2l 配准结果')

