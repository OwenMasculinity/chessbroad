clc;
clear all
close all

%% 边缘提取
Ref_im = rgb2gray(imread('Ref.png'));
Ref_im=double(Ref_im);
[M,N]=size(Ref_im);

figure;imshow(Ref_im,[]);title('Ref.png input');

sobelKernelY = [1 2 1;0 0 0; -1 -2 -1];
sobelKernelX = [-1 0 1;-2 0 2; -1 0 1];

derivativeX2 = imfilter(Ref_im, sobelKernelX,'replicate');
derivativeY2 = imfilter(Ref_im, sobelKernelY,'replicate');
gradientMagnitude_Ref = sqrt(derivativeX2.^2 + derivativeY2.^2);
figure;imshow(gradientMagnitude_Ref,[]);title('Ref.png边缘提取后结果');


%% 图像膨胀
gradientMagnitude_Ref = imdilate(gradientMagnitude_Ref,[0 1 0;1 1 1;0 1 0]);

%% 二值化
segmentationResultByMy=gradientMagnitude_Ref;
for i=1:M
    for j=1:N
        segmentationResultByMy(i,j)=(gradientMagnitude_Ref(i,j)>55);
    end
end

%% 中值滤波
segmentationResultByMy=medfilt2(segmentationResultByMy,[3,3]);
figure;imshow(segmentationResultByMy,[]);title('Ref.png二值化后结果');


%% hough变化:直线检测
% 获取格点间隔，用于final.m中参数设置
[H,T,R] = hough(segmentationResultByMy);
figure;
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');title('Ref.png图像hough变换后结果');
axis on, axis normal, hold on;

P  = houghpeaks(H,30,'threshold',ceil(0.3*max(H(:))));

x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');

lines = houghlines(segmentationResultByMy,T,R,P,'FillGap',10,'MinLength',100);
figure, imshow(segmentationResultByMy,[]), title('Ref.png直线提取'); hold on

i=1;
j=1;
for k = 1:length(lines)
    if(lines(k).theta==0)
        V_lines(i)=lines(k);
        i=i+1;
    end
    if(lines(k).theta==-90)
        H_lines(j)=lines(k);
        j=j+1;
    end
end

for i=1:length(V_lines)
    V_xy(:,:,i) = [V_lines(i).point1; V_lines(i).point2];
    plot(V_xy(:,1,i),V_xy(:,2,i),'LineWidth',2,'Color','green');
end

for i=1:length(H_lines)
    H_xy(:,:,i) = [H_lines(i).point1; H_lines(i).point2];
    plot(H_xy(:,1,i),H_xy(:,2,i),'LineWidth',2,'Color','red');
end
