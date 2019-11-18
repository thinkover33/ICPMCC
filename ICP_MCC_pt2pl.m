function [ R, t, corr, Error_of_iter, TData, TotalStep, TotalTime ] = ICP_MCC_pt2pl( DData, MData, R, t, MDataNormals)
%ICP_MCC_PT2PL 此处显示有关此函数的摘要
% This is an implementation of the Iterative Closest Point algorithm based on maximum correntropy criterion.
% The function takes two data sets(2D or 3D) and registers DData with MData.
% The code iterates till no more correspondences can be found.
% DData->MData
% Arguments: MData - M x D matrix of the x, y and z coordinates of data set 1, D is the dimension of point set
%            DData - N x D matrix of the x, y and z coordinates of data set 2
%            MDataNormals - The normals for each point in MData
%
% Returns: R - D x D accumulative rotation matrix used to register DData
%          t - D x 1 accumulative translation vector used to register DData
%          corr - M x D matrix of the index no.s of the corresponding points of
%                 MData and DData and their corresponding Euclidean distance
%          error - the  error between the corresponding points
%                  of MData and DData (normalized with res)
%          TData - M x D matrix of the registered DData
% 对于三维点集：
% 求解参数为 x = [α β γ tx ty tz]'
% 对于二维点集：
% 求解参数为 x = [α tx ty]'

% 参数设置
Epsilon = 10^-7;
MaxStep = 30                 ; % the max step of iteration
IS_Step_Show = 1                   ;
sigma_times = 80                ; %控制初始的sigma值是meanDist的多少倍
DecayPram = 0.97                   ; %控制sigma的值是否在循环中改变
isShowIni =1                  ;%控制是否显示初值
neighborNum = 10               ; %计算法向量选择的最近邻点数目

[nD, ~] = size(DData);   % the number of the model point
[nM, Dim] = size(MData);   % the number of the test point
if Dim<2
    error('Demension error! Only 2D or 3D data sets are supported.');
end
if Dim>3
    error('Demension error! Only 2D or 3D data sets are supported.');
end

if nargin <= 2
    R = eye(Dim);
end
if nargin <= 3
    TData = (R * DData')';
    mX = mean(TData);
    mY = mean(MData);
    t = (mY - mX)';
%       t = zeros(Dim,1);
end
if nargin <= 4
    if Dim ==3
        ptMData = pointCloud(MData);
        MDataNormals = pcnormals(ptMData,neighborNum);
    end
    if Dim == 2
        MDataNormals = PCANormals(MData);
    end
end


% Initialization
totaltime = 0;
TotalStep = 1;  % last step, if the program is pause, it is the total iterative step
CurrStep = 1;  % current step
LastMinuserror = 10^6; %Save the last and last step's error minus last step's error
Error_of_iter = zeros(MaxStep,1);
Error_of_iter(1) = 10^6;   % The Error
% 评估点集内点的密集程度：求其最近点的平均距离  评估结果为 meanDist
[~,TD] = knnsearch(DData,DData,'k', 2);
meanDist = mean(TD(:,2));
midDist = median(TD(:,2));
% 由点集的密度决定sigma的值
sigma0 = midDist * sigma_times;
sigma1 = sigma0;
% To obtain the transformation data
TData = (R * DData')'+ repmat(t',[nD,1]);
if isShowIni
    ICPPlot(TData,MData);
    title('初值设定');
end

% 函数回显进度
dotnum=1;
fprintf('ICP_MCC_pt2pl is processing')
if IS_Step_Show == 1
    gcf_step = figure;
end
if Dim == 3
    eulXYZ = rotm2eul(R);
    x = [eulXYZ t']' ; % 待求解量
    while (LastMinuserror > Epsilon && TotalStep < MaxStep)
        tic;
%         [corr, ~] = knnsearch(MData, TData, 'k',1);
        [corr, ~] = multiQueryKNNSearchImpl(pointCloud(MData), TData,1);
        for ii = 1:5
            A =  [cross(TData,MDataNormals(corr,:),2),MDataNormals(corr,:) ]; %此处到底是用TData 还是 DData 存疑
            b = dot(MDataNormals(corr,:) , (MData(corr,:)-TData),2);
            g = exp( -((A*x-b).^2)./ (2*sigma1^2) );
            x = (A.*repmat(g,[1,6]))'*A \ A' .* repmat(g',[6,1]) *b;
            
            cx = cos(x(1));
            cy = cos(x(2));
            cz = cos(x(3));
            sx = sin(x(1));
            sy = sin(x(2));
            sz = sin(x(3));
            R1 = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
                cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
                -sy,          sx*cy,          cx*cy];
            T1 = x(4:6);
            
            R = R1*R;
            t = R1*t + T1;
            TData = (R * DData')';
            TData = TData + repmat(t',[nD,1]);
        end
        tempdistance = (TData-MData(corr, :));
        Error_of_iter(CurrStep) = sum(sum(tempdistance.^2, 2))/nD;  % compute error
        
        TotalStep = CurrStep;     % TotalStep record current step as last step
        CurrStep = CurrStep + 1;  % CurrStep is next step
        % sigma的退火
        if sigma1 > meanDist * 2
            sigma1 = sigma1 * DecayPram;
        end
        
        if TotalStep == 1
            LastMinuserror = Error_of_iter(TotalStep);
        else
            LastMinuserror = abs(Error_of_iter(TotalStep) - Error_of_iter(TotalStep-1));
        end
        passtime = toc;
        % 函数回显
        if (TotalStep/10)>dotnum
            fprintf('.');
            dotnum = dotnum + 1;
        end
        %   disp(['Complete ', num2str(TotalStep),  ' time, ', 'take ', num2str(passtime), 's' ]);
        totaltime = totaltime + passtime;
        if IS_Step_Show == 1
            clf(gcf_step);
            ICPPlot(TData,MData,gcf_step);
            view(0,90);
            title(['Frame',num2str(CurrStep)]);
            %         waitforbuttonpress;
            pause(0.001);
        end
        
    end
else
    x = [acos(R(1,1)) t']';
    while(LastMinuserror > Epsilon && TotalStep < MaxStep)
        tic;
        [corr, ~] = knnsearch(MData, TData, 'k',1);
        detmat = TData(:,1).*MDataNormals(corr,2) - TData(:,2).*MDataNormals(corr,1);
        A = [detmat MDataNormals(corr,:)];
        b = dot(MDataNormals(corr,:) , (MData(corr,:)-TData),2);
        g = exp( -((A*x-b).^2)./ (2*sigma1^2) );
        x = (A.*repmat(g,[1,3]))'*A \ A' .* repmat(g',[3,1]) *b;
        
        R1 = [cos(x(1)) -sin(x(1)); sin(x(1)) cos(x(1))];
        T1 = x(2:3);
        R = R1*R;
        t = R1*t + T1;
        TData = (R * DData')';
        TData = TData + repmat(t',[nD,1]);
        tempdistance = (TData-MData(corr(:,1), :));
        Error_of_iter(CurrStep) = sum(sum(tempdistance.^2, 2))/nD;  % compute error
        
        TotalStep = CurrStep;     % TotalStep record current step as last step
        CurrStep = CurrStep + 1;  % CurrStep is next step
        % sigma的退火
        if sigma1 > meanDist * 2
            sigma1 = sigma1 * DecayPram;
        end
        
        if TotalStep == 1
            LastMinuserror = Error_of_iter(TotalStep);
        else
            LastMinuserror = abs(Error_of_iter(TotalStep) - Error_of_iter(TotalStep-1));
        end
        % 函数回显
        if (TotalStep/10)>dotnum
            fprintf('.');
            dotnum = dotnum + 1;
        end
        passtime = toc;
        %   disp(['Complete ', num2str(TotalStep),  ' time, ', 'take ', num2str(passtime), 's' ]);
        totaltime = totaltime + passtime;
        if IS_Step_Show == 1
            clf(gcf_step);
            ICPPlot(TData,MData,gcf_step);
            view(0,90);
            title(['Frame',num2str(CurrStep)]);
            %         waitforbuttonpress;
            pause(0.001);
        end
    end
end
if IS_Step_Show == 1
    close(gcf_step);
end
% 函数回显
fprintf('\n');
disp(['Complete ', num2str(TotalStep),  ' times.  Totaltime: ', num2str(totaltime), 'S.']);
disp(['meanDist: ', num2str(meanDist), ' midDist: ', num2str(midDist), '  Sigma: ' ,num2str(sigma0), ' -> ', num2str(sigma1)])
TotalTime = totaltime;
% disp(['Totaltime: ', num2str(totaltime), 'S' ]);
% disp(['Sigma:',num2str(sigma1)]);
end

