% Hemmorhages
I = imread('C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\Image processing\16_left.jpeg');
greenC = I(:,:,2);
comp = imcomplement(greenC);
histe = adapthisteq(comp);
adjustImage = imadjust(histe,[],[],6);
comp = imcomplement(adjustImage);
J = imadjust(comp,[],[],4);
J = imcomplement(J);
J = imadjust(J,[],[],9);
J = im2bw(J,0.2);
BW2 = bwareaopen(J, 7000);
wname = 'sym4';
[CA,CH,CV,CD] = dwt2(BW2,wname,'mode','per');
figure,imshow(CA),title('first image');
[lab num] = bwlabeln(CA, 8);
imtool(lab)
% % disp ('Area of Hemmorhages');
% areaa = bwarea(CA);
% disp(areaa);

b = bwboundaries(CA);
% axes(handles.axes5);
I = imresize(I,[500 752]);
figure,imshow(I)

hold on

% for k = 1:numel(b)
%     plot(b{k}(:,2), b{k}(:,1), 'w', 'Linewidth', 2)
% end
% imsave

