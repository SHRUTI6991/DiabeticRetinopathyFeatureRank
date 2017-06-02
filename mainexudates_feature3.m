%% Extract optic disc and artifacts from one image
tic
% Read image
retinaRGB = imread('C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\Image processing\21.jpg');
% Resize image
retinaRGB = resizeretina(retinaRGB, 752, 500);
% Get optic disc mask
closingThresholdValue = 0.64;
opticDiscDilationSize = 4;
artifactMinSize = 1100;
%[opticDiscMask, artifactsMask] = getopticdiscartifacts(retinaRGB, ...
           % closingThresholdValue, opticDiscDilationSize, artifactMinSize);
        
toc

%% Extract exudates from one image
tic
fileName = 'C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\Image processing\21.jpg';    
% Read image
retinaRGB = imread(fileName);
% Resize image
retinaRGB = resizeretina(retinaRGB, 752, 500);
% Read optic disc mask
opticDiscMask = imread(strrep(fileName, '.jpg', '.png'));
xx = double(opticDiscMask);
z = mesh(xx);
artifactsMask = imread(strrep(fileName, '.jpg', '.png'));
xx1 = double(artifactsMask);
z1 = mesh(xx1);
% Get optic disc mask
opticDiscDilation = 10;

%% Get intensity
    % subplot(1, 2, 1), imshow(retinaRGB); title('RGB');
    I = double(retinaRGB) / 255;
    I = sum(I, 3) ./ 3;
    % subplot(1, 2, 2), imshow(I); title('Intensity');

    %% Median filter on intensity channel
    % subplot(1, 2, 1), imshow(I); title('Before median filter');
    I = medfilt2(I);
    %I = imgaussfilt(I);
    % subplot(1, 2, 2), imshow(I); title('Median filter on intensity channel');

    %% Histogram equalization
    % subplot(1, 2, 1), imshow(I); title('Before histogram equalization');
    I = adapthisteq(I);
    % subplot(1, 2, 2), imshow(I); title('Histogram equalization');    

  
    %% Remove vessels by grayscale closing
    
    % subplot(1, 2, 1), imshow(I); title('Before grayscale closing');
    se = strel('disk', 8);
    closeI = imclose(I, se);
    
    %% Local standard deviation of an image
    % subplot(1, 2, 1), imshow(closeI); title('Before standard deviation');
    deviation = stdfilt(closeI, ones(7)); 
    % subplot(1, 2, 2), imshow(deviation, []); title('Standard deviation');

    %% Threshold and dilation
    % subplot(1, 2, 1), imshow(deviation, []); title('Before threshold and dilation');
    level = graythresh(deviation);
    mask = im2bw(deviation, level);
    se = strel('disk', 6);
    mask = imdilate(mask, se);
    % subplot(1, 2, 2), imshow(mask); title('Thresholded and dilated');
    
      %% Watershed Transformation
    
    %mask = watershed(mask);
    

    %% Create region of interest
    retinaMask = im2bw(I, 0.2);
    retinaMask = imfill(retinaMask, 'holes');
    se = strel('disk', 16);
    retinaMask = imerode(retinaMask, se);
    % subplot(1, 2, 1), imshow(I), title('Original image');
    % subplot(1, 2, 2), imshow(I, 'InitialMag', 'fit')
    % Make a truecolor all-green image.
    % green = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
    % hold on
    % h = imshow(green);
    % hold off
    % Use our influence map as the AlphaData for the solid green image.
    % set(h, 'AlphaData', retinaMask)

    %% Remove circular shape around retina
    % subplot(1, 2, 1), imshow(mask); title('Before removing circular shape around');
    maskOfCenter = mask .* retinaMask;
    % subplot(1, 2, 2), imshow(maskOfCenter); title('Only region of interest');
    
    %% Flood fill
    % subplot(1, 2, 1), imshow(maskOfCenter); title('Before filling');
    maskFilled = imfill(maskOfCenter, 'holes');
    size(maskFilled);
    % subplot(1, 2, 2), imshow(maskFilled); title('Flood filled');

    %% Remove optic disc
    % subplot(1, 2, 1), imshow(maskFilled); title('Before optic disc elimination');
    se = strel('disk', opticDiscDilation);
    opticDiscMask = imdilate(z, se);
    maskOfInterest = uint8(retinaRGB) .* imcomplement(z);
    
     %% Remove artifacts
    % subplot(1, 2, 1), imshow(maskOfInterest); title('Before artifacts elimination');
    se = strel('disk', opticDiscDilation);
    artifactsMask = imdilate(z1, se);
    maskOfInterest = uint8(maskOfInterest) .* imcomplement(z1);
    xx2 = double(maskOfInterest);
    z3 = mesh(xx2);
    % subplot(1, 2, 2), imshow(maskOfInterest); title('Without artifacts');
    
    %% Overlay mask on the original image
    % subplot(1, 2, 1), imshow(I); title('Before overlay');
    marker = I .* imcomplement(z3);
    % subplot(1, 2, 2), imshow(marker); title('Overlay');

    %% Reconstruction
    % subplot(1, 2, 1), imshow(marker); title('Before reconstruction');
    reconstructed = imreconstruct(marker, I);
    % subplot(1, 2, 2), imshow(reconstructed); title('Reconstruction');

    %% Threshold on image differences
    diff = I - reconstructed;
    % subplot(1, 2, 1), imshow(diff, []), title('Difference before threshold');
    % level = graythresh(diff)
    level = 0.01;
    exudatesMask = im2bw(diff, level);
    % subplot(1, 2, 2), imshow(exudatesMask), title('Mask');

    %% Overlay exudates mask on the original image
    subplot(1, 2, 1), imshow(retinaRGB), title('Original image');
    subplot(1, 2, 2), imshow(reconstructed, 'InitialMag', 'fit')
%     % Make a truecolor all-green image.
%     green1 = cat(3, zeros(size(I)), ones(size(I)), zeros(size(I)));
%     
%     hold on
%     h = imshow(green1);
%     hold off
%     % Use our influence map as the AlphaData for the solid green image.
%     set(h, 'AlphaData', exudatesMask)



%exudatesMask = getexudates(retinaRGB, opticDiscMask, artifactsMask, opticDiscDilation);
toc

% %% Postprocess one image
% tic
% fileName = 'C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\DRIVE\DRIVE\training\images\21_training.tif';  % 225_left
% % Read images
% retina = imread(fileName);
% exudates = imread(strrep(fileName, '.jpeg', '_exudates.png'));
% % Postprocessing
% exudatesMaxSize = 160;
% exudatesPostprocessed = postprocessing(exudates, retina, exudatesMaxSize);
% 
% toc

% %% Extract features from one image
% % Load red lesions, exudates and optic disc mask
% fileName = 'C:\Users\SHRUTI\Desktop\Project\diabetic-retinopathy-master\DRIVE\DRIVE\training\images\21_training.tif';  % 225_left
% exudates = imread(strrep(fileName, '.jpeg', '_exudates.png'));
% opticDisc = imread(strrep(fileName, '.jpeg', '_optic_disc_mask.png'));
% redLesions = imread(strrep(strrep(fileName, '.jpeg', '_redlesion.png'), ...
%                             'E:/Dev/CAD/Diabetic Retinopathy/train', ...
%                             'E:/Dropbox/Diabetic Retinopathy/red_lesions_all'));
% % Make the images logical
% opticDisc = im2bw(opticDisc, 0.1);
% redLesions = im2bw(redLesions, 0.1);
% 
% % Exudates postprocessing
% retina = imread(fileName);
% exudatesMaxSize = 160;
% exudates = postprocessing(exudates, retina, exudatesMaxSize);
% 
% % Features 1 - 12 for exudates
% featuresExudates = getlesionsfeatures(exudates, opticDisc)
% % Features 1 - 12 (13 - 24) for red lesions
% featuresRedLesions = getlesionsfeatures(redLesions, opticDisc)
% % 12. Optic disc distance from center
% opticDistance = getopticdistance(opticDisc)
% 
% %% Extract all optic discs and artifacts and save them as files
% % parpool(2);  % Set workers for paraller computations
% % tic
% 
% % List all images
% directory = 'E:/Dev/CAD/Diabetic Retinopathy/test/';
% filesLeft = dir(strcat(directory, '*_left.jpeg'));
% filesRight = dir(strcat(directory, '*_right.jpeg'));
% files = [filesLeft; filesRight];
% % For each image
% nFiles = length(files);
% for i = 30000  : nFiles
%     fileName = strcat(directory, files(i).name);
%     fprintf('Optic disc and artifacts, processing image %i / %i, %s.\n', i, nFiles, fileName);
%     
%     % Read image
%     retinaRGB = imread(fileName);
%     % Resize image
%     retinaRGB = resizeretina(retinaRGB, 752, 500);
%     % Get optic disc mask
%     closingThresholdValue = 0.64;
%     opticDiscDilationSize = 4;
%     artifactMinSize = 1100;
%     [opticDiscMask, artifactsMask] = getopticdiscartifacts(retinaRGB, ...
%             closingThresholdValue, opticDiscDilationSize, artifactMinSize);
%     
%     % Save the mask
%     imwrite(opticDiscMask, strrep(fileName, '.jpeg', '_optic_disc_mask.png'))
%     imwrite(artifactsMask, strrep(fileName, '.jpeg', '_artifacts_mask.png'))
%     
% end
% % toc
% % delete(gcp);  % Close threads pool
% 
% %% Extract all exudates and save them as files
% % parpool(2);  % Set workers for paraller computations
% % tic
% 
% % List all images
% directory = 'E:/Dev/CAD/Diabetic Retinopathy/test/';
% filesLeft = dir(strcat(directory, '*_left.jpeg'));
% filesRight = dir(strcat(directory, '*_right.jpeg'));
% files = [filesLeft; filesRight];
% % For each image
% nFiles = length(files);
% for i = 25000 : nFiles
%     fileName = strcat(directory, files(i).name);
%     fprintf('Exudates, processing image %i / %i, %s.\n', i, nFiles, fileName);
%     
%     % Read image
%     retinaRGB = imread(fileName);
%     % Resize image
%     retinaRGB = resizeretina(retinaRGB, 752, 500);
%     % Read optic disc mask
%     opticDiscMask = imread(strrep(fileName, '.jpeg', '_optic_disc_mask.png'));
%     artifactsMask = imread(strrep(fileName, '.jpeg', '_artifacts_mask.png'));
%     % Get optic disc mask
%     opticDiscDilation = 10;
%     exudatesMask = getexudates(retinaRGB, opticDiscMask, artifactsMask, opticDiscDilation);
%     
%     % Save the exudates
%     imwrite(exudatesMask, strrep(fileName, '.jpeg', '_exudates.png'))
%     
% end
% % toc
% % delete(gcp);  % Close threads pool
% 
% %% Extract all features from train data and save them to file
% parpool(2);  % Set workers for paraller computations
% tic
% 
% % Get all file names
% names = readtable('D:/CAD/diabetic-retinopathy/DRClassification/trainLabels.csv', ...
%                   'Delimiter', ',', 'ReadVariableNames', true);
% names = names.image;
% % Initialize features matrix
% nImages = size(names, 1);
% nFeatures = 33;
% features = zeros(nImages, nFeatures);
% 
% % For each image
% for i = 1 : nImages
%     if (mod(i, 100) == 0)
%         fprintf('Features train, processing image %i / %i, %s.\n', i, nImages, names{i});
%     end
%     
%     % Read images
%     exudates = imread(strcat('E:/Dev/CAD/Diabetic Retinopathy/train/', names{i}, '_exudates.png'));
%     opticDisc = imread(strcat('E:/Dev/CAD/Diabetic Retinopathy/train/', names{i}, '_optic_disc_mask.png'));
%     redLesions = imread(strcat('E:/Dropbox/Diabetic Retinopathy/red_lesions_all/', names{i}, '_redlesion.png'));
%     
%     % Make the images logicalp
%     opticDisc = im2bw(opticDisc, 0.1);
%     redLesions = im2bw(redLesions, 0.1);
%     
%     % Exudates postprocessing
%     % retina = imread(strcat('E:/Dev/CAD/Diabetic Retinopathy/train/', names{i}, '.jpeg'));
%     % exudatesMaxSize = 160;
%     % exudates = postprocessing(exudates, retina, exudatesMaxSize);
%     
%     % Features 1 - 11 for exudates
%     featuresExudates = getlesionsfeatures(exudates, opticDisc);
%     % Features 1 - 11 (12 - 22) for red lesions
%     featuresRedLesions = getlesionsfeatures(redLesions, opticDisc);
%     % 12. Optic disc distance from center
%     opticDistance = getopticdistance(opticDisc);
%     
%     % Save features into matrix
%     features(i,:) = [featuresRedLesions, featuresExudates, opticDistance];
%     
% end
% 
% % Write csv file
% csvwrite('D:/CAD/diabetic-retinopathy/DRClassification/features_all_postprocessed.csv', features);
% 
% toc
% delete(gcp);  % Close threads pool
% 
% %% Extract all features from test data and save them to file
% % parpool(2);  % Set workers for paraller computations
% tic
% 
% % List all images
% directory = 'E:/Dev/CAD/Diabetic Retinopathy/test/';
% filesLeft = dir(strcat(directory, '*_left.jpeg'));
% filesRight = dir(strcat(directory, '*_right.jpeg'));
% files = [filesLeft; filesRight];
% 
% % Get all file names
% % names = readtable('D:/CAD/diabetic-retinopathy/DRClassification/trainLabels.csv', ...
% %                   'Delimiter', ',', 'ReadVariableNames', true);
% % names = names.image;
% 
% % Initialize features matrix
% nImages = numel(files);
% nFeatures = 33;
% features = zeros(nImages, nFeatures);
% names = cell(nImages, 1);
% 
% % For each image
% for i = 1 : nImages
%     fileName = strcat(directory, files(i).name);
%     fileNameShort = strrep(files(i).name, '.jpeg', '');
%     
%     if (mod(i, 100) == 0)
%         fprintf('Features test, processing image %i / %i, %s.\n', i, nImages, fileName);
%     end
%     
%     % Read images
%     exudates = imread(strcat('E:/Dev/CAD/Diabetic Retinopathy/test/', fileNameShort, '_exudates.png'));
%     opticDisc = imread(strcat('E:/Dev/CAD/Diabetic Retinopathy/test/', fileNameShort, '_optic_disc_mask.png'));
%     redLesions = imread(strcat('E:/Dev/CAD/Diabetic Retinopathy/red_lesions_test/', fileNameShort, '_redlesion.png'));
%     
%     % Make the images logicalp
%     opticDisc = im2bw(opticDisc, 0.1);
%     redLesions = im2bw(redLesions, 0.1);
%     
%     % Features 1 - 11 for exudates
%     featuresExudates = getlesionsfeatures(exudates, opticDisc);
%     % Features 1 - 11 (12 - 22) for red lesions
%     featuresRedLesions = getlesionsfeatures(redLesions, opticDisc);
%     % 12. Optic disc distance from center
%     opticDistance = getopticdistance(opticDisc);
%     
%     % Save features into matrix and names intro array
%     features(i,:) = [featuresRedLesions, featuresExudates, opticDistance];
%     names{i} = fileNameShort;
%     
% end
% 
% % Write csv file with features
% csvwrite('D:/CAD/diabetic-retinopathy/DRClassification/features_test.csv', features);
% 
% % Write csv file with names
% fid = fopen('D:/CAD/diabetic-retinopathy/DRClassification/names.csv', 'w') ;
% for i = 1 : nImages
%     fprintf(fid, '%s\n', names{i, end});
% end
% fclose(fid) ;
% 
% toc
% % delete(gcp);  % Close threads pool
