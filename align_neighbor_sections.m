clearvars -except I
close all
addpath ./tiff_loading/utilities
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar');
myimfuse = @(x,y)(imfuse(x,y,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));

sliceA=3;
sliceB=4;
slide=1;
gcamp_channel=2;
path='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah025/';
files=dir(path);

t=0;
for i=1:length(files)
    tmp=strsplit(files(i).name,'.');
    if strcmpi(tmp{end},'nd2');
        t=t+1;
        filename{t}=[files(t).folder '/' files(i).name];
    end
end

if ~exist('I');
    for i=1:length(filename)
        I{i}=double(load_tiff(filename{i}));
        I{i}=I{i}(:,:,:,gcamp_channel);
        info=nd2finfo(filename{i});
        tmp=zeros(round(info.img_height*info.calib_factor),round(info.img_width*info.calib_factor),size(I{i},3));
        for z=1:size(I{i},3)
            tmp(:,:,z)=imresize(I{i}(:,:,z),[round(info.img_height*info.calib_factor) round(info.img_width*info.calib_factor)]);
            [i z]
        end
        I{i}=tmp;
    end
end

dims=zeros(1,3);
for i=1:length(I)
    dims=max([dims; size(I{i})],[],1);
end

for i=1:length(I)
    I{i}=padarray(I{i},[dims(1)-size(I{i},1) dims(2)-size(I{i},2) 0],'post');
end



I1slice=max(I{sliceA}(:,:,:),[],3);
I2slice=max(I{sliceB}(:,:,:),[],3);%
I2slice_warped=I2slice;

notGood=1;
rigid=1;
affine=0;
globTform=[eye(2) zeros(2,1);zeros(1,2) 1];
while notGood==1
    close all
    imagesc(10*imfuse(I1slice,I2slice_warped,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
    title('Click successive pairs of matching landmarks - Chose RED first then GREEN  - right click when done');
    
    [x,y]=getpts;
    
    x(end)=[];
    y(end)=[];
    
    matches1=[x(mod(1:length(x),2)==1) y(mod(1:length(x),2)==1)];
    matches2=[x(mod(1:length(x),2)==0) y(mod(1:length(x),2)==0)];
    
    %     subplot(1,2,1)
    %     imagesc(I1slice*1000);
    %     hold on
    %     for i=1:size(matches1,1)
    %         plot(matches1(i,1),matches1(i,2),'r.','MarkerSize',20);
    %     end
    %
    %     subplot(1,2,2)
    %     imagesc(I2slice*1000);
    %     hold on
    %     for i=1:size(matches2,1)
    %         plot(matches2(i,1),matches2(i,2),'r.','MarkerSize',20);
    %     end
    
    if rigid==1
        [R,T]=wahba(matches2,matches1);
        tform=affine2d([R zeros(2,1);T 1]);
    elseif affine==1
        beta=linsolve([matches2 ones(size(matches2,1),1)],matches1);
        tform=affine2d([beta(1:2,:) zeros(2,1);beta(3,:) 1]);
    end
    
    globTform=tform.T*globTform;
    
    % tform = cp2tform(detections2(matching_cells(:,2),:),detections1(matching_cells(:,1),:),'affine');
    
    
    
    I2slice_warped=imwarp(I2slice,affine2d(globTform),"OutputView",imref2d(size(I1slice)));
    close all
    imagesc(10*imfuse(I1slice,I2slice_warped,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
    switch input('Is it good enough? (1 for yes) Do you want to switch to affine? (2 for affine)? Do you want to switch to rigid? (3 for rigid)');
        case 1
            notGood=0;
            
        case 2
            rigid=0;
            affine=1;
        case 3
            rigid=1;
            affine=0;
    end
    
    
end

radius=40;
sigma=200;

close all
tform=affine2d(globTform);
I2slice_globtform=imwarp(I2slice,tform,"OutputView",imref2d(size(I1slice)));
imagesc(10*myimfuse(I2slice_globtform,I1slice));

for i=1:size(I{sliceB},3)
    I2transformed(:,:,i)=imwarp(I{sliceB}(:,:,i),tform,"OutputView",imref2d(size(I1slice)));
end
I1=I{sliceA};
I2=I{sliceB};
detections_moving=myDetections(I2transformed);
detections_fixed=myDetections(I1);

[idx,D]=knnsearch(detections_fixed,detections_moving);


for i=1:size(detections_moving,1)
    if D(i)<radius
        hold on
        plot(detections_moving(i,1),detections_moving(i,2),'y.','MarkerSize',10);
        plot(detections_fixed(idx(i),1),detections_fixed(idx(i),2),'c.','MarkerSize',10);
        
        plot([detections_moving(i,1) detections_fixed(idx(i),1)],[detections_moving(i,2) detections_fixed(idx(i),2)],'w-','LineWidth',2);
    end
end


% vfieldTotal=pointDeformer(I2transformed(:,:,8),I1(:,:,1),sigma,detections_moving(D<radius,:),detections_fixed(idx(D<radius),:));
vfieldTotal=pointDeformerIterative(I2transformed(:,:,8),I1(:,:,1),sigma);
close all

for i=1:size(I{sliceB},3)
    I2transformed_warped(:,:,i)=imwarp(I2transformed(:,:,i),vfieldTotal);
end


[ax, ~] = tight_subplot(1, 3, [0.01 0.01], 0.1,0.1);
% ax1=subplot(1,3,1);
axes(ax(1));
imagesc(10*myimfuse(I2(:,:,8),I1(:,:,1)));title('Unregistered');axis equal;axis off;text(100,100,'Unregistered','color','w','FontWeight','bold','FontSize',20);
% ax2=subplot(1,3,2);
axes(ax(2))
imagesc(10*myimfuse(I2transformed(:,:,8),I1(:,:,1)));title('Rigid/Affine');axis equal;axis off;text(100,100,'Rigid/affine','color','w','FontWeight','bold','FontSize',20);
% ax3=subplot(1,3,3);
axes(ax(3))
imagesc(10*myimfuse(I2transformed_warped(:,:,8),I1(:,:,1)));title('Deformable');axis equal;axis off;text(100,100,'Deformable','color','w','FontWeight','bold','FontSize',20);
linkaxes([ax(1) ax(2) ax(3)]);
set(gcf,'color','w');


% ax1=subplot(1,3,1);
% imagesc(10*myimfuse(I2(:,:,8),I1(:,:,1)));title('Unregistered');axis off;
% ax2=subplot(1,3,2);
% imagesc(10*myimfuse(I2transformed(:,:,8),I1(:,:,1)));title('Rigid/Affine');axis off
% ax3=subplot(1,3,3);
% imagesc(10*myimfuse(I2transformed_warped(:,:,8),I1(:,:,1)));title('Deformable');axis off
% linkaxes([ax1 ax2 ax3]);
% set(gcf,'color','w');
% save(['slide_' num2str(slide) '_slice_' num2str(sliceB) '_to_slice_' num2str(sliceA) '_deformable_tform.mat'],'globTform','vfieldTotal');



