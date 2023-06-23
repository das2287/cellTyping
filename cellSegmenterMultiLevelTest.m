function Xtest=cellSegmenterMultiLevelTest(I,levels)
for ds=1:length(levels)
    for z=1:size(I,3)
        I_DS{ds}(:,:,z)=imresize(imresize(I(:,:,z),levels(ds),'bilinear'),[size(I,1) size(I,2)],'bilinear');
    end
end


shifts=[1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
Xtest=zeros(numel(I_DS{1}),size(shifts,1)*size(I_DS,2));
for d=1:size(I_DS,2)
    for s=1:size(shifts,1)
        Ishift=imtranslate(I_DS{d},shifts(s,:));
        Xtest(:,(d-1)*size(shifts,1)+s)=Ishift(:);
    end
end
end