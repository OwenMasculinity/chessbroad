clc;
clear all;
close all;


%% ��Ե��ȡ
Test_im = rgb2gray(imread('Test.png'));
im1=double(Test_im);
[M, N] = size(Test_im);

Ref_im = rgb2gray(imread('Ref.png'));
im2=double(Ref_im);

figure;imshow(im1,[]);title('Test.png input');
figure;imshow(im2,[]);title('Ref.png input');

sobelKernelY = [1 2 1;0 0 0; -1 -2 -1];
sobelKernelX = [-1 0 1;-2 0 2; -1 0 1];

derivativeX1 = imfilter(im1, sobelKernelX,'replicate');
derivativeY1 = imfilter(im1, sobelKernelY,'replicate');
gradientMagnitude_Test = sqrt(derivativeX1.^2 + derivativeY1.^2);

derivativeX2 = imfilter(im2, sobelKernelX,'replicate');
derivativeY2 = imfilter(im2, sobelKernelY,'replicate');
gradientMagnitude_Ref = sqrt(derivativeX2.^2 + derivativeY2.^2);

figure;imshow(gradientMagnitude_Test,[]);title('Test.png��Ե��ȡ�����');
% dot_Test=ginput();       %ȡ�ĸ��㣬���������ϣ����ϣ����£�����
figure;imshow(gradientMagnitude_Ref,[]);title('Ref.png��Ե��ȡ����');
% dot_Ref=ginput();

%% ȡ�ĸ����̽ǵ����꣬����Corner_xy.m�ó��������ǣ����ϣ����ϣ����£����£�
dot_Test=[761.5,103.5;1725.5,343.5;141.5,769.5;1441.5,1225.5];%Test.png�ĸ����̽ǵ�����
dot_Ref=[94.5,438.5;1238.5,442.5;102.5,1562.5;1238.5,1562.5]; %Ref.png�ĸ����̽ǵ�����

%% ͸��任
%�ĸ�ԭ����
y=[dot_Test(1,1) dot_Test(2,1) dot_Test(3,1) dot_Test(4,1)];       
x=[dot_Test(1,2) dot_Test(2,2) dot_Test(3,2) dot_Test(4,2)];

%�ĸ��¶���
Y=[dot_Ref(1,1) dot_Ref(2,1) dot_Ref(3,1) dot_Ref(4,1)];     
X=[dot_Ref(1,2) dot_Ref(2,2) dot_Ref(3,2) dot_Ref(4,2)];

B=[X(1) Y(1) X(2) Y(2) X(3) Y(3) X(4) Y(4)]';%�任����ĸ����㣬�����ұߵ�ֵ
%�����ⷽ���飬���̵�ϵ��
A=[x(1) y(1) 1 0 0 0 -X(1)*x(1) -X(1)*y(1);             
0 0 0 x(1) y(1) 1 -Y(1)*x(1) -Y(1)*y(1);
   x(2) y(2) 1 0 0 0 -X(2)*x(2) -X(2)*y(2);
0 0 0 x(2) y(2) 1 -Y(2)*x(2) -Y(2)*y(2);
   x(3) y(3) 1 0 0 0 -X(3)*x(3) -X(3)*y(3);
0 0 0 x(3) y(3) 1 -Y(3)*x(3) -Y(3)*y(3);
   x(4) y(4) 1 0 0 0 -X(4)*x(4) -X(4)*y(4);
0 0 0 x(4) y(4) 1 -Y(4)*x(4) -Y(4)*y(4)];

fa=inv(A)*B;        %���ĵ���õķ��̵Ľ⣬Ҳ��ȫ�ֱ任ϵ��
a=fa(1);b=fa(2);c=fa(3);
d=fa(4);e=fa(5);f=fa(6);
g=fa(7);h=fa(8);

rot=[d e f;
     a b c;
     g h 1];        %��ʽ�е�һ������x,Matlab��һ����ʾy�������Ҿ���1,2�л�����

pix1=rot*[1 1 1]'/(g*1+h*1+1);  %�任��ͼ�����ϵ�
pix2=rot*[1 N 1]'/(g*1+h*N+1);  %�任��ͼ�����ϵ�
pix3=rot*[M 1 1]'/(g*M+h*1+1);  %�任��ͼ�����µ�
pix4=rot*[M N 1]'/(g*M+h*N+1);  %�任��ͼ�����µ�

height=round(max([pix1(1) pix2(1) pix3(1) pix4(1)])-min([pix1(1) pix2(1) pix3(1) pix4(1)]));     %�任��ͼ��ĸ߶�
width=round(max([pix1(2) pix2(2) pix3(2) pix4(2)])-min([pix1(2) pix2(2) pix3(2) pix4(2)]));      %�任��ͼ��Ŀ��
imgn=zeros(height,width);

delta_y=round(abs(min([pix1(1) pix2(1) pix3(1) pix4(1)])));            %ȡ��y����ĸ��ᳬ����ƫ����
delta_x=round(abs(min([pix1(2) pix2(2) pix3(2) pix4(2)])));            %ȡ��x����ĸ��ᳬ����ƫ����
inv_rot=inv(rot);

for i = 1-delta_y:height-delta_y                        %�ӱ任ͼ���з���Ѱ��ԭͼ��ĵ㣬������ֿն�������ת�Ŵ�ԭ��һ��
    for j = 1-delta_x:width-delta_x
        pix=inv_rot*[i j 1]';       %��ԭͼ�������꣬��Ϊ[YW XW W]=fa*[y x 1],�������������[YW XW W],W=gy+hx+1;
        pix=inv([g*pix(1)-1 h*pix(1);g*pix(2) h*pix(2)-1])*[-pix(1) -pix(2)]'; %�൱�ڽ�[pix(1)*(gy+hx+1) pix(2)*(gy+hx+1)]=[y x],����һ�����̣���y��x�����pix=[y x];
        
        if pix(1)>=0.5 && pix(2)>=0.5 && pix(1)<=M && pix(2)<=N
            imgn(i+delta_y,j+delta_x)=gradientMagnitude_Test(round(pix(1)),round(pix(2)));     %���ڽ���ֵ,Ҳ������˫���Ի�˫������ֵ
        end  
    end
end

imgn=imgn'; %����ͼ�η���
figure;
imshow(uint8(imgn));title('Test.png͸��仯����');
figure;imhist(uint8(imgn));title('Test.png͸��任���ֱ��ͼ');

%% ͼ������
im_temp = uint8(imgn);
im_temp = imdilate(im_temp,[0 1 0;1 1 1;0 1 0]);
imshow(im_temp);title('Test.pngͼ�����ͺ���');

%% ��ֵ��
[M,N]=size(im_temp);
segmentationResultByMy=im_temp;
for i=1:M
    for j=1:N
        segmentationResultByMy(i,j)=(im_temp(i,j)>55);
    end
end

segmentationResultByMy=medfilt2(segmentationResultByMy,[1,1]);%��ֵ�˲�
figure;imshow(segmentationResultByMy,[]);title('Test.pngͼ���ֵ������');


%% hough�仯:ֱ�߼��
[H,T,R] = hough(segmentationResultByMy);
figure;
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');title('Test.pngͼ��hough�任����');
axis on, axis normal, hold on;

P  = houghpeaks(H,24,'threshold',ceil(0.3*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');

% ֱ����ȡ
lines = houghlines(segmentationResultByMy,T,R,P,'FillGap',10,'MinLength',100);
figure, imshow(segmentationResultByMy,[]), title('ֱ����Բ��ȡ');hold on
max_len = 0;

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

%% ȷ�����̽���
%�ο�Ref.pngͼƬ���������ȡ��1��1������λ�ã������̸����
min_lines=114;%�������������̺���֮�����ĳ�ֵ
min_columns=125;%���������������̺���֮�����ĳ�ֵ
point=zeros(10,9,2);

%��1,1��λ��Ϊ��1014,1035��
% ��10*9�����
for i=1:10
    for j=1:9
        point(i,j,:)=[1014+(i-1)*min_lines,1035+(j-1)*min_columns];
        plot(point(i,j,1),point(i,j,2),'x','LineWidth',2,'Color','blue');
    end
end

% %% ���Բ��
% [centers, radii] = imfindcircles(im_temp,[45 60],'ObjectPolarity','dark', ...
%     'Sensitivity',0.95,'Method','twostage','EdgeThreshold',0.14);
% [m,n]=size(centers);
% 
% 
% %% �������λ�ò���ʾ����
% disp('����λ��:');
% for i=1:10
%     for j=1:9
%         for p=1:m
%             if((centers(p,1)<point(i,j,1)+min_lines/4)&&(centers(p,1)>point(i,j,1)-min_lines/2.5))
%                 if((centers(p,2)<point(i,j,2)+min_columns/4)&&(centers(p,2)>point(i,j,2)-min_columns/2.5))
%                     %h = viscircles(centers(p,:),radii(p,:));
%                     s=sprintf('(%d,%d)',i,j);  
%                     disp(s);
%                     break;
%                 end
%             end
%         end
%     end
% end
% 
% 
% %% ��ɫȷ��
% im_temp_r = im2bw(imread('red.png'));
% im_temp_snapshot = imread('gray2bw.png');
% figure;imshow(im_temp_r)
% disp('����ɫ��Ϣ�������ַ�λ����Ϣ���:');
% temp_black = zeros(100,2);
% temp_red = zeros(100,2);
% Temp = zeros(100,4);
% snapshot = zeros(71,71,32);
% n = 1;
% k = 1;
% l = 1;
% for i=1:10
%     for j=1:9
%         for p=1:m
%             if((centers(p,1)<point(i,j,1)+min_lines/4)&&(centers(p,1)>point(i,j,1)-min_lines/2.5))
%                 if((centers(p,2)<point(i,j,2)+min_columns/4)&&(centers(p,2)>point(i,j,2)-min_columns/2.5))
%                     h = viscircles(centers(p,:),radii(p,:));
%                     r = round(centers(p,1));
%                     c = round(centers(p,2));
%                     snapshot(:,:,n) = imcrop(im_temp_snapshot,[r-35, c-35, 70, 70]);
%                     for t = -10 : 10
%                         Temp(n,1) = Temp(n,1) + double(im_temp_r(c,r+t));
%                     end
%                     if Temp(n,1) < 15
%                         s = sprintf('Black(%d,%d)',i,j);
%                         temp_black(k,:) = [r, c];
%                         k = k + 1;
%                         Temp(n,2) = ('B');
%                         Temp(n,3:4) = [i, j];
%                     else
%                         s = sprintf('Red(%d,%d)',i,j);
%                         temp_red(l,:) = [r, c];
%                         l = l + 1;
%                         Temp(n,2) = ('R');
%                         Temp(n,3:4) = [i, j];
%                     end
%                     n = n + 1;
% %                     disp(s);
%                     break;
%                 end
%             end
%         end
%     end
% end
% 
% %% �����ַ�ȷ����ģ��ƥ��
% % Snapshot characters in circles above
% templatePath='./templateCharacter/';
% fileFormat='.bmp';
% templateImage=zeros(20,20,32);
% Timage=zeros(32,800);
% for i=1:32 
%     stri=num2str(i);
%     imagePath=[templatePath,stri,fileFormat];
%     tempImage=im2bw(imread(imagePath));
%     templateImage(:,:,i)=tempImage;
%     clear imagePath stri tempImage;
% end
% characterImage=zeros(20,20,32);
% for i = 1:32
%     characterImage(:,:,i) = im2bw(imresize(snapshot(:,:,i),[20,20]));
% end
% characterImage = 1 - characterImage;
% Uimage=zeros(32,800);
% 
% 
% Y=zeros(1,32);
% for i=1:32
%     U=length(find( characterImage(:,:,i))~=0);
%     for j=1:32
%         T=length(find( templateImage(:,:,j))~=0);
%         tempV=characterImage(:,:,i)& templateImage(:,:,j);
%         V=length(find(tempV)~=0);
%         tempW=xor(tempV,templateImage(:,:,j));
%         W=length(find(tempW)~=0);
%         tempX=xor(tempV,characterImage(:,:,i));
%         X=length(find(tempX)~=0);
%         TUV=(T+U+V)/3;
%         tempSum=sqrt(((T-TUV)*(T-TUV)+(U-TUV)*(U-TUV)+(V-TUV)*(V-TUV))/2);
%         Y(j)=V/(W/T*X/U*tempSum);
%     end
%     [MAX,indexMax]=max(Y);
%     stri=num2str(indexMax);
%     switch stri
%         case '1'
% %             disp(strcat(Temp(i,2),'Ju',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ju(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '2'
% %             disp(strcat(Temp(i,2),'Ma',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ma(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '3'
% %             disp(strcat(Temp(i,2),'XiangB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,XiangB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '4'
% %             disp(strcat(Temp(i,2),'ShiB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,ShiB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '5'
% %             disp(strcat(Temp(i,2),'Jiang',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Jiang(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '6'
% %             disp(strcat(Temp(i,2),'ShiB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,ShiB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '7'
% %             disp(strcat(Temp(i,2),'XiangR',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,XiangR(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '8'
% %             disp(strcat(Temp(i,2),'Ma',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ma(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '9'
% %             disp(strcat(Temp(i,2),'Ju',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ju(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '10'
% %             disp(strcat(Temp(i,2),'PaoB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,PaoB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '11'
% %             disp(strcat(Temp(i,2),'PaoR',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,PaoR(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '12'
% %             disp(strcat(Temp(i,2),'Zu',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Zu(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '13'
% %             disp(strcat(Temp(i,2),'Bing',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Bing(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '14'
% %             disp(strcat(Temp(i,2),'Zu',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Zu(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '15'
% %             disp(strcat(Temp(i,2),'Bing',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Bing(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '16'
% %             disp(strcat(Temp(i,2),'Zu',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Zu(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '17'
% %             disp(strcat(Temp(i,2),'Bing',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Bing(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '18'
% %             disp(strcat(Temp(i,2),'Zu',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Zu(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '19'
% %             disp(strcat(Temp(i,2),'Bing',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Bing(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '20'
% %             disp(strcat(Temp(i,2),'Zu',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Zu(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '21'
% %             disp(strcat(Temp(i,2),'Bing',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Bing(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '22'
% %             disp(strcat(Temp(i,2),'PaoR',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,PaoR(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '23'
% %             disp(strcat(Temp(i,2),'PaoB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,PaoB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '24'
% %             disp(strcat(Temp(i,2),'Ju',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ju(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '25'
% %             disp(strcat(Temp(i,2),'Ma',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ma(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '26'
% %             disp(strcat(Temp(i,2),'XiangR',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,XiangR(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '27'
% %             disp(strcat(Temp(i,2),'ShiB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,ShiB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '28'
% %             disp(strcat(Temp(i,2),'Shuai',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Shuai(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '29'
% %             disp(strcat(Temp(i,2),'ShiR',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,ShiR(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '30'
% %             disp(strcat(Temp(i,2),'XiangB',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,XiangB(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '31'
% %             disp(strcat(Temp(i,2),'Ma',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ma(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%         case '32'
% %             disp(strcat(Temp(i,2),'Ju',Temp(i,3),Temp(i,4)));
%             str = sprintf('%s,Ju(%d,%d)',Temp(i,2),Temp(i,3),Temp(i,4));
%     end
%     disp(str);
%     clear imagePath indexMax;
% end
