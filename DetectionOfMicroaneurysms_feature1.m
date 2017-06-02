I = imread('C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\Image processing\16_left.jpeg');
	
greenC = I(:,:,2);
comp = imcomplement(greenC);
histe = adapthisteq(comp);
adjustImage = imadjust(histe,[],[],3);
comp = imcomplement(adjustImage);
J = imadjust(comp,[],[],4);
J = imcomplement(J);
J = imadjust(J,[],[],4);
K=fspecial('disk',5);
L=imfilter(J,K,'replicate');
L = im2bw(L,0.4);
M =  bwmorph(L,'tophat');
% M = im2bw(M);
wname = 'sym4';
[CA,CH,CV,CD] = dwt2(M,wname,'mode','per');
figure,imshow(CA);

b = bwboundaries(CA);
I = imresize(I,[303 350]);
%figure, imshow(I);
hold on

% corners = detectFASTFeatures(I,'MinContrast',0.1);
% J = insertMarker(I,corners,'circle');


for area_bloodvessels = 1:numel(b)
    plot(b{area_bloodvessels}(:,2), b{area_bloodvessels}(:,1), 'b', 'Linewidth', 1)
end 
