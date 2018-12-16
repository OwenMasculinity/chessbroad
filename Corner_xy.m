clc;
clear all;
close all;

%% ��ͼ
im_Ref=imread('Ref.png');
im_Ref=rgb2gray(im_Ref);
imshow(im_Ref);title('Ref.png input');
[M_Ref,N_Ref]=size(im_Ref);
im_Ref=double(im_Ref);

im_Test=imread('Test.png');
im_Test=rgb2gray(im_Test);
figure;imshow(im_Test);title('Test.png input');
[M_Test,N_Test]=size(im_Test);
im_Test=double(im_Test);

%% ��Ե��ȡ
sobelKernelY = [1 2 1;0 0 0; -1 -2 -1];
sobelKernelX = [-1 0 1;-2 0 2; -1 0 1];

derivativeX2_Ref = imfilter(im_Ref, sobelKernelX,'replicate');
derivativeY2_Ref = imfilter(im_Ref, sobelKernelY,'replicate');
im_Ref = sqrt(derivativeX2_Ref.^2 + derivativeY2_Ref.^2);

derivativeX2_Test = imfilter(im_Test, sobelKernelX,'replicate');
derivativeY2_Test = imfilter(im_Test, sobelKernelY,'replicate');
im_Test = sqrt(derivativeX2_Test.^2 + derivativeY2_Test.^2);
% figure;imshow(im_Ref,[]);
% figure;imshow(im_Test,[]);

%% ��ֵ��
for i=1:M_Ref
    for j=1:N_Ref
        im_Ref(i,j)=(im_Ref(i,j)>55);
    end
end

for i=1:M_Test
    for j=1:N_Test
        im_Test(i,j)=(im_Test(i,j)>55);
    end
end

%% ��ֵ�˲�
im_Ref=medfilt2(im_Ref,[1,1]);
im_Test=medfilt2(im_Test,[1,1]);


%% ��ȡ�ǵ�
% for Ref.png
[H_Ref,T_Ref,R_Ref] = hough(im_Ref);
P_Ref  = houghpeaks(H_Ref,30,'threshold',ceil(0.3*max(H_Ref(:))));
lines_Ref = houghlines(im_Ref,T_Ref,R_Ref,P_Ref,'FillGap',15,'MinLength',500);

a=10000;
b=0;
for k = 1:length(lines_Ref)
   if(lines_Ref(k).theta==0)
       if(lines_Ref(k).rho<a)
          a=lines_Ref(k).rho;
          XY1_Ref=[lines_Ref(k).point1; lines_Ref(k).point2];
       end
       if(lines_Ref(k).rho>b)
          b=lines_Ref(k).rho;
          XY2_Ref=[lines_Ref(k).point1; lines_Ref(k).point2];
       end
   end
end

XY1_1=[XY1_Ref(1,1),XY1_Ref(1,2)];%����
XY1_2=[XY1_Ref(2,1),XY1_Ref(2,2)];%����
XY1_3=[XY2_Ref(2,1),XY2_Ref(2,2)];%����
XY1_4=[XY2_Ref(2,2),XY2_Ref(2,2)];%����

figure, imshow(im_Ref,[]), hold on
plot(XY1_Ref(1,1),XY1_Ref(1,2),'o','LineWidth',2,'Color','red');
plot(XY1_Ref(2,1),XY1_Ref(2,2),'o','LineWidth',2,'Color','red');
plot(XY2_Ref(1,1),XY2_Ref(1,2),'o','LineWidth',2,'Color','red');
plot(XY2_Ref(2,1),XY2_Ref(2,2),'o','LineWidth',2,'Color','red');

% ͼƬRef.png ���̽ǵ�����
XY_Ref=[XY1_1;XY1_3;XY1_2;XY1_4]; %˳��Ϊ�������ϣ����ϣ����£����£�


% for Test.png
[H_Test,T_Test,R_Test] = hough(im_Test);
P_Test  = houghpeaks(H_Test,30,'threshold',ceil(0.3*max(H_Test(:))));
lines_Test = houghlines(im_Test,T_Test,R_Test,P_Test,'FillGap',15,'MinLength',500);

a=10000;
b=0;
for k = 1:length(lines_Test)
   if(lines_Test(k).theta>0)
       if(lines_Test(k).rho<a)
          a=lines_Test(k).rho;
          XY1_Test=[lines_Test(k).point1; lines_Test(k).point2];
       end
       if(lines_Test(k).rho>b)
          b=lines_Test(k).rho;
          XY2_Test=[lines_Test(k).point1; lines_Test(k).point2];
       end
   end
end

XY2_1=[XY1_Test(1,1),XY1_Test(1,2)];%����
XY2_2=[XY1_Test(2,1),XY1_Test(2,2)];%����
XY2_3=[XY2_Test(2,1),XY2_Test(2,2)];%����
XY2_4=[XY2_Test(2,2),XY2_Test(2,2)];%����

figure, imshow(im_Test,[]), hold on
plot(XY1_Test(1,1),XY1_Test(1,2),'o','LineWidth',2,'Color','red');
plot(XY1_Test(2,1),XY1_Test(2,2),'o','LineWidth',2,'Color','red');
plot(XY2_Test(1,1),XY2_Test(1,2),'o','LineWidth',2,'Color','red');
plot(XY2_Test(2,1),XY2_Test(2,2),'o','LineWidth',2,'Color','red');

% ͼƬ Test.png ���̽ǵ�����
XY_Test=[XY2_1;XY2_3;XY2_2;XY2_4];%˳��Ϊ�������ϣ����ϣ����£����£�

%% ��ʾ���̽ǵ�����
disp('Ref.png���̽ǵ����꣺');
disp(XY_Ref);
disp('Test.png���̽ǵ����꣺');
disp(XY_Test);