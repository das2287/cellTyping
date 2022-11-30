clearvars -except I1 Iz1 props1 Iz_t_labels1 im1 imz1 imz_t1 roiImage1 I2 Iz2 props2 Iz_t_labels2 im2 imz2 imz_t2 roiImage2
clc
close all
vec=@(x)(x(:));
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
addpath ./tiff_loading/utilities
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar');


% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice1_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice1_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice2_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice2_Round2_PV_SATB1002.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice3_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice3_Round2_PV_SATB1.nd2';

round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice4_Round1.nd2';
round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice4_Round2_PV_SATB1.nd2';


round_1_shortname=strsplit(round_1_filename,'/');round_1_shortname=round_1_shortname{end};
round_2_shortname=strsplit(round_2_filename,'/');round_2_shortname=round_2_shortname{end};

max_sigma=50;%%Maximum size of cells
min_sigma=30;%%Minimum size of cells
z_threshold=2;%%Cell brightness threshold
volumeThreshold=300;%%Cell size threshold
patchRadius=25; %% size of patches to show cells at the last step
gcamp_channel=2;%Which channel is Gcamp
% colors=[0 0.5 1;0 1 0;1 0 0;0 1 1];
colors=[0 1 1;1 0 0;0 1 0;0 0.5 1];%% color map of 4 channels
%% Step 0 - load datasets
if ~exist('I1')
    I1=double(load_tiff(round_1_filename));
end

if ~exist('I2')
    I2=double(load_tiff(round_2_filename));
end

round_1_dims=input('Enter pixel dimensions for Round 1 [x y z]','s');
round_1_dims=str2num(round_1_dims);

round_2_dims=input('Enter pixel dimensions for Round 2 [x y z]','s');
round_2_dims=str2num(round_2_dims);


if ~exist('im1')
im1=[repmat(minmax(max(I1(:,:,:,1),[],3)),[1 1 3]).*reshape(colors(1,:),[1 1 3]) repmat(minmax(max(I1(:,:,:,2),[],3)),[1 1 3]).*reshape(colors(2,:),[1 1 3]);
    repmat(minmax(max(I1(:,:,:,3),[],3)),[1 1 3]).*reshape(colors(3,:),[1 1 3]) repmat(minmax(max(I1(:,:,:,4),[],3)),[1 1 3]).*reshape(colors(4,:),[1 1 3])];
end

if ~exist('im2')
im2=[repmat(minmax(max(I2(:,:,:,1),[],3)),[1 1 3]).*reshape(colors(1,:),[1 1 3]) repmat(minmax(max(I2(:,:,:,2),[],3)),[1 1 3]).*reshape(colors(2,:),[1 1 3]);
    repmat(minmax(max(I2(:,:,:,3),[],3)),[1 1 3]).*reshape(colors(3,:),[1 1 3]) repmat(minmax(max(I2(:,:,:,4),[],3)),[1 1 3]).*reshape(colors(4,:),[1 1 3])];
end
% figure(1);
% imagesc(im*5);
% drawnow

if ~exist('roiImage1')
    roiImage1=minmax(squeeze(max(I1(:,:,:,2:4),[],3)))*25;
end

if ~exist('roiImage2')
    roiImage2=minmax(squeeze(max(I2(:,:,:,2:4),[],3)))*25;
end

figure

imagesc([roiImage1 roiImage2]);
title('Draw ROI in left image for registerable cells - then repeat for right image');
h1 = drawfreehand;
h2 = drawfreehand;

hold on
plot(h1.Position(:,1),h1.Position(:,2),'g*');
plot(h2.Position(:,1),h2.Position(:,2),'g*');

h2.Position(:,1)=h2.Position(:,1)-size(roiImage1,1);

im1_roi=h1.Position;
im2_roi=h2.Position;


figure
subplot(1,2,1)
imagesc(roiImage1)
hold on
plot(im1_roi(:,1),im1_roi(:,2),'g*');

subplot(1,2,2)
imagesc(roiImage2)
hold on
plot(im2_roi(:,1),im2_roi(:,2),'g*');






%% Step 1 - Detect cells

if ~exist('Iz1')
    Iz1=zeros(size(I1,1),size(I1,2),size(I1,3));
    for ch=gcamp_channel
        Iz1=spatial_zscore(I1(:,:,:,ch),max_sigma)-spatial_zscore(I1(:,:,:,ch),min_sigma);    
    end
end

if ~exist('Iz2')
    Iz2=zeros(size(I2,1),size(I2,2),size(I2,3));
    for ch=gcamp_channel
        Iz2=spatial_zscore(I2(:,:,:,ch),max_sigma)-spatial_zscore(I2(:,:,:,ch),min_sigma);    
    end
end


Iz_t1=Iz1.*(Iz1>z_threshold);
Iz_t2=Iz2.*(Iz2>z_threshold);

if ~exist('Iz_t_labels1');
    tic;Iz_t_labels1=bwlabeln(Iz_t1,26);toc
end

if ~exist('Iz_t_labels2');
    tic;Iz_t_labels2=bwlabeln(Iz_t2,26);toc
end

if ~exist('props1')
    tic;props1=regionprops3(Iz_t_labels1);toc
end

if ~exist('props2')
    tic;props2=regionprops3(Iz_t_labels2);toc
end


    loc1=props1.Centroid;
    vol1=props1.Volume;
    
    loc2=props2.Centroid;
    vol2=props2.Volume;
    
    %% cross register here
    
    %% transform image and centroids here
    
    
    in1=inpolygon(loc1(:,1),loc1(:,2),im1_roi(:,1),im1_roi(:,2));
    in2=inpolygon(loc2(:,1),loc2(:,2),im2_roi(:,1),im2_roi(:,2));
    detections1=loc1(and(in1==1,vol1>volumeThreshold),:);
    detections2=loc2(and(in2==1,vol2>volumeThreshold),:);

figure;
subplot(1,2,1)
detectionsImage1=zeros(size(I1,1),size(I1,2),3);
detectionsImage1(:,:,1)=minmax(max(I1(:,:,:,gcamp_channel),[],3));

imagesc(detectionsImage1*5);
hold on
for i=1:size(detections1,1)
    text(detections1(i,1),detections1(i,2),num2str(i),'Color','w','FontWeight','bold','FontSize',15);
end
plot(detections1(:,1),detections1(:,2),'wo','MarkerSize',10);
title('Round 1');

subplot(1,2,2)
detectionsImage2=zeros(size(I2,1),size(I2,2),3);
detectionsImage2(:,:,1)=minmax(max(I2(:,:,:,gcamp_channel),[],3));

imagesc(detectionsImage2*5);
hold on
for i=1:size(detections2,1)
    text(detections2(i,1),detections2(i,2),num2str(i),'Color','w','FontWeight','bold','FontSize',15);
end
plot(detections2(:,1),detections2(:,2),'wo','MarkerSize',10);
title('Round 2');
drawnow


matching_cells=input('Enter matching cells in Round 1 and Round 2 images (input format ([1 1;2 2;3 3])');


beta=linsolve(detections2(matching_cells(:,2),1:2),detections1(matching_cells(:,1),1:2));
beta=blkdiag(beta,1);
tform=affine3d([beta zeros(3,1);zeros(1,3) 1]);

I2warped=imwarp(I2(:,:,:,gcamp_channel),tform,'outputview',imref3d(size(I1(:,:,:,gcamp_channel))));


figure('units','normalized','outerposition',[0 0 1 1/2])
subplot(1,3,1);
imagesc(detectionsImage1*5);
title(round_1_shortname,'Interpreter','none');
subplot(1,3,2)
imagesc(5*imfuse(max(I2warped,[],3),max(I1(:,:,:,gcamp_channel),[],3)));
title('Cross registration result')
subplot(1,3,3);
imagesc(detectionsImage2*5);
title(round_2_shortname,'Interpreter','none');
set(gcf,'color','w');



[optimizer,metric] = imregconfig("multimodal");
moving_reg = imregister(max(I2warped,[],3),max(I1(:,:,:,gcamp_channel),[],3),'affine',optimizer,metric);

figure('units','normalized','outerposition',[0 0 1 1/2])
subplot(1,3,1);
imagesc(detectionsImage1*5);
title(round_1_shortname,'Interpreter','none');
subplot(1,3,2)
imagesc(imfuse(max(moving_reg,[],3),max(I1(:,:,:,gcamp_channel),[],3)));
title('Cross registration result')
subplot(1,3,3);
imagesc(detectionsImage2*5);
title(round_2_shortname,'Interpreter','none');
set(gcf,'color','w');

% %% Step 2 - cell type
% t=0;
% for ch=gcamp_channel
%     for i=1:size(detections{ch},1)
%         t=t+1;
%         detectionsConsolidated(t,:)=detections{ch}(i,:);
%     end
% end
% 
% 
% 
% t=0;
% for i=1:size(detectionsConsolidated,1)
%     for ch=1:4
%         try;
%             tmp=max(I1(detectionsConsolidated(i,2)-patchRadius:detectionsConsolidated(i,2)+patchRadius,detectionsConsolidated(i,1)-patchRadius:detectionsConsolidated(i,1)+patchRadius,:,ch),[],3);
%             if ch==1
%                 t=t+1;
%             end
%             patch{t,ch}=tmp;
%         end
%     end
% end
% 
% num_per_row=floor(sqrt(size(patch,1)))*4;
% 
% remaining=num_per_row-mod(size(patch,1),num_per_row);
% if mod(remaining,num_per_row)==0
%     remaining=0;
% end
% 
% sizePatch=size(patch,1);
% for t=1:remaining
%     for ch=1:4
%     patch{sizePatch+t,ch}=ones(size(patch{1,ch}));
%     end
% end
% 
% num_in_row=0;
% stackPicture=[];
% barcodePicture=[];
% for i=1:size(patch,1)
%     
%     tmp=[];
%     for ch=1:4
%         tmp=[tmp minmax(repmat(patch{i,ch},[1 1 3])).*reshape(colors(ch,:),[1 1 3])];
%     end
%     barcodePicture=[barcodePicture;tmp];
%     num_in_row=num_in_row+1;
%     if num_in_row==num_per_row
%         num_in_row=0;
%         stackPicture=[stackPicture ones(size(barcodePicture,1),5,3) barcodePicture];
%         barcodePicture=[];
%     end
% end
% 
% figure
% imagesc(permute(stackPicture,[2 1 3]));axis equal;axis off
% set(gcf,'color','w')
% title(['ROI Volume: ' num2str(polyarea(h.Position(:,1)*round_1_dims(1),h.Position(:,2)*round_1_dims(2))*round_1_dims(3)) ' um^3']);
%% Step 3 - write to excel



