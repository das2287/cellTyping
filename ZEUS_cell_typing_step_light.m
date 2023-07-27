clear all
clc
close all
vec=@(x)(x(:));
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));
myimfuse = @(x,y)(imfuse(x,y,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
addpath ./tiff_loading/utilities
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar');

% 
round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice1_Round1.nd2';
round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice1_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice2_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice2_Round2_PV_SATB1002.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice3_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice3_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice4_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice4_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice5_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice5_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice6_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice6_Round2_PV_SATB1.nd2';


round_1_shortname=strsplit(round_1_filename,'/');round_1_shortname=round_1_shortname{end};
round_2_shortname=strsplit(round_2_filename,'/');round_2_shortname=round_2_shortname{end};

max_sigma=50;%%Maximum size of cells
min_sigma=30;%%Minimum size of cells
z_threshold=0.5;%%Cell brightness threshold
volumeThreshold=50;%%Cell size threshold
patchRadius=25; %% size of patches to show cells at the last step
gcamp_channel=2;%Which channel is Gcamp
colors=[0 1 1;1 0 0;0 1 0;0 0.5 1;0 1 1;1 0 0;0 1 0;0 0.5 1];%% color map of 4 channels
dims=[1 1]; %X-Y pixel dimensions
%% Step 0 - load datasets
if ~exist('I0')
    I0{1}=double(load_tiff(round_1_filename));
    I0{2}=double(load_tiff(round_2_filename));
end

if ~exist('I0warped')
    a=strsplit(round_2_filename,'.');
    load([a{1} '_deformation_map.mat'],'vfieldTotal','vfield_final','globTform');
    for ch=1:size(I0{2},4)
        tmp=I0{2}(:,:,:,ch);
        for z=1:size(I0{2},3)
            I0warped{2}(:,:,z,ch)=imwarp(imwarp(tmp(:,:,z),affine2d(globTform),"OutputView",imref2d(size(I0{1}(:,:,z,1)))),vfieldTotal+vfield_final);
        end
    end
end
if ~exist('I')
    
    I=cat(4,max(I0{1}(:,:,1:min(size(I0{1},3),size(I0{2},3)),:),[],3),max(I0warped{2}(:,:,1:min(size(I0{1},3),size(I0{2},3)),:),[],3));
end


if ~exist('im')
    im=[repmat(minmax(max(I(:,:,:,1),[],3)),[1 1 3]).*reshape(colors(1,:),[1 1 3]) repmat(minmax(max(I(:,:,:,2),[],3)),[1 1 3]).*reshape(colors(2,:),[1 1 3]);
        repmat(minmax(max(I(:,:,:,3),[],3)),[1 1 3]).*reshape(colors(3,:),[1 1 3]) repmat(minmax(max(I(:,:,:,4),[],3)),[1 1 3]).*reshape(colors(4,:),[1 1 3])];
end
% figure(1);
% imagesc(im*5);
% drawnow


%% Step 1 - Detect cells

if ~exist('Iz')
    Iz=zeros(size(I));
    for ch=1:size(I,4)
        Iz(:,:,:,ch)=imgaussfilt(spatial_zscore(I(:,:,:,ch),max_sigma)-spatial_zscore(I(:,:,:,ch),min_sigma),1);
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
        tic;props{ch}=regionprops('table',Iz_t_labels(:,:,:,ch));toc
    end
end


for ch=1:4
    loc=props{ch}.Centroid;
    vol=props{ch}.Area;
    detections{ch}=loc(vol>volumeThreshold,:);
    area{ch}=vol(vol>volumeThreshold);
end



%% Step 2 - cell type
clear detectionsConsolidated patch patch_location
t=0;
for ch=gcamp_channel
    for i=1:size(detections{ch},1)
        t=t+1;
        detectionsConsolidated(t,:)=detections{ch}(i,:);
    end
end



t=0;
detectionsRemoved=zeros(size(detectionsConsolidated,1),1);
for i=1:size(detectionsConsolidated,1)
    for ch=1:8
        try;
            tmp=max(I(detectionsConsolidated(i,2)-patchRadius:detectionsConsolidated(i,2)+patchRadius,detectionsConsolidated(i,1)-patchRadius:detectionsConsolidated(i,1)+patchRadius,:,ch),[],3);
            if ch==1
                t=t+1
            end
            patch{t,ch}=tmp;patch_location{t}=detectionsConsolidated(i,:);
        catch
            detectionsRemoved(i)=1;
        end
    end
end

detectionsConsolidated(detectionsRemoved==1,:)=[];


num_per_row=floor(sqrt(size(detectionsConsolidated,1)))*4;

remaining=num_per_row-mod(size(patch,1),num_per_row);
if mod(remaining,num_per_row)==0
    remaining=0;
end

sizePatch=size(detectionsConsolidated,1);
for t=1:remaining
    for ch=1:8
        patch{sizePatch+t,ch}=ones(size(patch{1,ch}));patch_location{sizePatch+t}=[nan nan];
    end
end

num_in_row=0;
stackPicture=[];
barcodePicture=[];
for i=1:size(patch,1)
    
    tmp=[];
    for ch=1:8
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

figure;
imagesc(10*myimfuse(max(I(:,:,:,gcamp_channel),[],3),max(I(:,:,:,gcamp_channel+4),[],3)));
hold on
for i=1:size(patch,1)
    text(patch_location{i}(1),patch_location{i}(2)-5,num2str(i),'Color','m','FontWeight','bold','FontSize',15);
end
title('Round 1');
drawnow





figure
imagesc(permute(stackPicture,[2 1 3]));axis equal;axis off
set(gcf,'color','w')
%% Step 3 - write to excel



