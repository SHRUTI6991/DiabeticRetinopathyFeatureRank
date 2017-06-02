%Program for Retinal Blood Vessel Extraction

%Author : Athi Narayanan S
%M.E, Embedded Systems,
%K.S.R College of Engineering
%Erode, Tamil Nadu, India.
%http://sites.google.com/site/athisnarayanan/
%s_athi1983@yahoo.co.in

%Program Description
%This program is the main entry of the application.
%This program extracts blood vessels from a retina image using Kirsch's Templates.

%Clear Memory & Command Window
clc;
clear all;
close all;

%Read Input Retina Image
inImg = imread('C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\Image processing\16_left.jpeg');
dim = ndims(inImg);
if(dim == 3)
    %Input is a color image
    inImg = rgb2gray(inImg);
end

%Extract Blood Vessels
Threshold = 10;
bloodVessels = VesselExtract(inImg, Threshold);

%Output Blood Vessels image

figure;
subplot(121);imshow(inImg);title('Input Image');
subplot(122);imshow(bloodVessels);title('Extracted Blood Vessels');