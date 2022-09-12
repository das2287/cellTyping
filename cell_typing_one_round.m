clearvars -except I Iz props Iz_t_labels im imz imz_t roiImage
clc
close all
vec=@(x)(x(:));
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
addpath ./tiff_loading/utilities
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar');

filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah102/sah102_section3_CCK_Cal_M2R_2x2.nd2';

max_sigma=50;%%Maximum size of cells
min_sigma=30;%%Minimum size of cells
z_threshold=2;%%Cell brightness threshold
volumeThreshold=300;%%Cell size threshold
patchRadius=25; %% size of patches to show cells at the last step

% colors=[0 0.5 1;0 1 0;1 0 0;0 1 1];
colors=[0 1 1;1 0 0;0 1 0;0 0.5 1];%% color map of 4 channels
%% Step 0 - load datasets
if ~exist('I')
    I=double(load_tiff(filename));
end

dims=input('Enter pixel dimensions [x y z]','s');
dims=str2num(dims);
if ~exist('im')
im=[repmat(minmax(max(I(:,:,:,1),[],3)),[1 1 3]).*reshape(colors(1,:),[1 1 3]) repmat(minmax(max(I(:,:,:,2),[],3)),[1 1 3]).*reshape(colors(2,:),[1 1 3]);
    repmat(minmax(max(I(:,:,:,3),[],3)),[1 1 3]).*reshape(colors(3,:),[1 1 3]) repmat(minmax(max(I(:,:,:,4),[],3)),[1 1 3]).*reshape(colors(4,:),[1 1 3])];
end
% figure(1);
% imagesc(im*5);
% drawnow

if ~exist('roiImage')
    roiImage=minmax(squeeze(max(I(:,:,:,2:4),[],3)))*25;
end
figure


imagesc(roiImage);
title('Draw ROI');
h = drawfreehand;
hold on
plot(h.Position(:,1),h.Position(:,2),'g*');

%% Step 1 - Detect cells

if ~exist('Iz')
    Iz=zeros(size(I));
    for ch=1:size(I,4)
        Iz(:,:,:,ch)=spatial_zscore(I(:,:,:,ch),max_sigma)-spatial_zscore(I(:,:,:,ch),min_sigma);    
    end
end


Iz_t=Iz.*(Iz>z_threshold);

if ~exist('Iz_t_labels');
for ch=1:4
    tic;Iz_t_labels(:,:,:,ch)=bwlabeln(Iz_t(:,:,:,ch),26);toc
end
end

if ~exist('props')
for ch=1:4
    tic;props{ch}=regionprops3(Iz_t_labels(:,:,:,ch));toc
end
end


for ch=1:4
    loc=props{ch}.Centroid;
    vol=props{ch}.Volume;
    in=inpolygon(loc(:,1),loc(:,2),h.Position(:,1),h.Position(:,2));
    detections{ch}=loc(and(in==1,vol>volumeThreshold),:);
end

figure;
imagesc(im*5);
hold on
plot(detections{1}(:,1),detections{1}(:,2),'wo','MarkerSize',10);
plot(detections{2}(:,1)+size(I,1),detections{2}(:,2),'wo','MarkerSize',10);
plot(detections{3}(:,1),detections{3}(:,2)+size(I,2),'wo','MarkerSize',10);
plot(detections{4}(:,1)+size(I,1),detections{4}(:,2)+size(I,2),'wo','MarkerSize',10);
title('Round 1');
drawnow


%% Step 2 - cell type
t=0;
for ch=3
    for i=1:size(detections{ch},1)
        t=t+1;
        detectionsConsolidated(t,:)=detections{ch}(i,:);
    end
end



t=0;
for i=1:size(detectionsConsolidated,1)
    for ch=1:4
        try;
            tmp=max(I(detectionsConsolidated(i,2)-patchRadius:detectionsConsolidated(i,2)+patchRadius,detectionsConsolidated(i,1)-patchRadius:detectionsConsolidated(i,1)+patchRadius,:,ch),[],3);
            if ch==1
                t=t+1;
            end
            patch{t,ch}=tmp;
        end
    end
end

num_per_row=floor(sqrt(size(patch,1)))*4;

remaining=num_per_row-mod(size(patch,1),num_per_row);
if mod(remaining,num_per_row)==0
    remaining=0;
end

sizePatch=size(patch,1);
for t=1:remaining
    for ch=1:4
    patch{sizePatch+t,ch}=ones(size(patch{1,ch}));
    end
end

num_in_row=0;
stackPicture=[];
barcodePicture=[];
for i=1:size(patch,1)
    
    tmp=[];
    for ch=1:4
        tmp=[tmp minmax(repmat(patch{i,ch},[1 1 3])).*reshape(colors(ch,:),[1 1 3])];
    end
    barcodePicture=[barcodePicture;tmp];
    num_in_row=num_in_row+1;
    if num_in_row==num_per_row
        num_in_row=0;
        stackPicture=[stackPicture ones(size(barcodePicture,1),5,3) barcodePicture];
        barcodePicture=[];
    end
end

figure
imagesc(permute(stackPicture,[2 1 3]));axis equal;axis off
set(gcf,'color','w')
title(['ROI Volume: ' num2str(polyarea(h.Position(:,1)*dims(1),h.Position(:,2)*dims(2))*dims(3)) ' um^3']);
%% Step 3 - write to excel



