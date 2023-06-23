clearvars -except I E
close all
addpath ./tiff_loading/utilities
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar');
myimfuse = @(x,y)(imfuse(x,y,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
minmax = @(x)((x-nanmin(x(:)))./(nanmax(x(:))-nanmin(x(:))));
vec=@(x)(x(:));
sliceA=3;
sliceB=4;
slide=1;
gcamp_channel=2;
path='/Users/erdem/Dropbox/FAST_datasets/immuno_exps/sah034/';
files=dir(path);

t=0;
for i=1:length(files)
    tmp=strsplit(files(i).name,'.');
    if strcmpi(tmp{end},'nd2');
        t=t+1;
        ex_vivo_filename{t}=[files(t).folder '/' files(i).name];
    end
end

if ~exist('I');
    for i=1:length(ex_vivo_filename)
        I{i}=double(load_tiff(ex_vivo_filename{i}));
        I{i}=I{i}(:,:,:,gcamp_channel);
        info=nd2finfo(ex_vivo_filename{i});
        tmp=zeros(round(info.img_height*info.calib_factor),round(info.img_width*info.calib_factor),size(I{i},3));
        for z=1:size(I{i},3)
            tmp(:,:,z)=imresize(I{i}(:,:,z),[round(info.img_height*info.calib_factor) round(info.img_width*info.calib_factor)],'bilinear');
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
    Ipad{i}=padarray(I{i},[dims(1)-size(I{i},1) dims(2)-size(I{i},2) 0],'post');
end

t=0;
for i=1:length(files)
    tmp=strsplit(files(i).name,'.');
    if strcmpi(tmp{end},'tiff');
        t=t+1;
        in_vivo_filename{t}=[files(t).folder '/' files(i).name];
    end
end

if ~exist('E');
    for i=1:length(in_vivo_filename)
        E{i}=double(load_tiff(in_vivo_filename{i}));
        info=imfinfo(in_vivo_filename{i});
        tmp=zeros(round(size(E{i},1)*info(1).XResolution),round(size(E{i},2)*info(1).YResolution),size(E{i},3));
        for z=1:size(E{i},3)
            tmp(:,:,z)=imresize(E{i}(:,:,z),[round(size(E{i},1)*info(1).XResolution) round(size(E{i},2)*info(1).YResolution)],'bilinear');
            [i z]
        end
        E{i}=tmp;
    end
end