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
% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice1_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice1_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice2_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice2_Round2_PV_SATB1002.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice3_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice3_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice4_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice4_Round2_PV_SATB1.nd2';

round_1_filename='/Users/admin/Dropbox/sah020_Slice1_Round1.nd2'
round_2_filename='/Users/admin/Dropbox/sah020_Slice1_Round2_PV_SATB1.nd2';

% round_1_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice6_Round1.nd2';
% round_2_filename='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah020/sah020_Slice6_Round2_PV_SATB1.nd2';


round_1_shortname=strsplit(round_1_filename,'/');round_1_shortname=round_1_shortname{end};
round_2_shortname=strsplit(round_2_filename,'/');round_2_shortname=round_2_shortname{end};

gcamp_channel=2;%Which channel is Gcamp
colors=[0 1 1;1 0 0;0 1 0;0 0.5 1];
%% Step 0 - load datasets
if ~exist('I0')
    I0{1}=double(load_tiff(round_1_filename));
    I0{2}=double(load_tiff(round_2_filename));
end
for i=1:length(I0)
    I{i}=I0{i}(:,:,:,gcamp_channel);
end

sliceA=1;
sliceB=2;

I1slice=max(I{sliceA},[],3);
I2slice=max(I{sliceB},[],3);
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
    
    if rigid==1
        [R,T]=wahba(matches2,matches1);
        tform=affine2d([R zeros(2,1);T 1]);
    elseif affine==1
        beta=linsolve([matches2 ones(size(matches2,1),1)],matches1);
        tform=affine2d([beta(1:2,:) zeros(2,1);beta(3,:) 1]);
    end
    
    globTform=tform.T*globTform;
    
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
sigma=50;

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


close all
[final_pass,vfield_final]=deformableReg(max(I1,[],3),max(I2transformed,[],3),20,50,200);

for i=1:size(I{sliceB},3)
    I2transformed_warped(:,:,i)=imwarp(I2transformed(:,:,i),vfield_final);
end


vfieldTotal=pointDeformerIterative(max(I2transformed_warped,[],3),max(I1,[],3),sigma);
close all

for i=1:size(I{sliceB},3)
    I2transformed_warped_final(:,:,i)=imwarp(I2transformed(:,:,i),vfieldTotal+vfield_final);
end



figure
[ax, ~] = tight_subplot(1, 4, [0.01 0.01], 0.1,0.1);
axes(ax(1));
imagesc(10*myimfuse(max(I2,[],3),max(I1,[],3)));title('Unregistered');axis equal;axis off;text(100,100,'Unregistered','color','w','FontWeight','bold','FontSize',20);
axes(ax(2))
imagesc(10*myimfuse(max(I2transformed,[],3),max(I1,[],3)));title('Rigid/Affine');axis equal;axis off;text(100,100,'Rigid/affine','color','w','FontWeight','bold','FontSize',20);
axes(ax(3))
imagesc(10*myimfuse(max(I2transformed_warped,[],3),max(I1,[],3)));title('Deformable Manual');axis equal;axis off;text(100,100,'Deformable Manual','color','w','FontWeight','bold','FontSize',20);

axes(ax(4))
imagesc(10*myimfuse(max(I2transformed_warped_final,[],3),max(I1,[],3)));title('Deformable Auto');axis equal;axis off;text(100,100,'Deformable Auto','color','w','FontWeight','bold','FontSize',20);
linkaxes([ax(1) ax(2) ax(3) ax(4)]);
set(gcf,'color','w');



for ch=1:size(I0{2},4)
    tmp=zscore(I0{2}(:,:,:,ch),[],'all');
    for z=1:size(I0{2},3)
        I0warped{2}(:,:,z,ch)=imwarp(imwarp(tmp(:,:,z),affine2d(globTform),"OutputView",imref2d(size(I1slice))),vfieldTotal+vfield_final);
    end
end

for ch=1:size(I0{1},4)
    tmp=zscore(I0{1}(:,:,:,ch),[],'all');
    for z=1:size(I0{1},3)
        I0warped{1}(:,:,z,ch)=tmp(:,:,z);
    end
end

figure
[ax, ~] = tight_subplot(2, 4, [0.01 0.01], 0.01,0.01);
for ch=1:4
    axes(ax(ch))
    imagesc(0.2*im2color(max(I0warped{1}(:,:,:,ch),[],3),colors(ch,:)));axis equal;axis off;
end

for ch=1:4
    axes(ax(4+ch))
    imagesc(0.2*im2color(max(I0warped{2}(:,:,:,ch),[],3),colors(ch,:)));axis equal;axis off;
end
linkaxes(ax(1:8));



a=strsplit(round_2_filename,'.');
save([a{1} '_deformation_map.mat'],'vfieldTotal','vfield_final','globTform');



