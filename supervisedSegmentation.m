function [Ipred,gpmdl]=supervisedSegmentation(I,levels)

for ds=1:length(levels)
    for z=1:size(I,3)
        I_DS{ds}(:,:,z)=imresize(imresize(I(:,:,z),levels(ds),'bilinear'),[size(I,1) size(I,2)],'bilinear');
        [ds z]
    end
end

[Xtrain,Ytrain,Xtest]=cellSegmenterMultiLevelTrain_v2(I,levels,100);
% Xtest=cellSegmenterMultiLevelTest(I,levels);

gpmdl=fitrgp(Xtrain,Ytrain);

Ypred=predict(gpmdl,Xtest);

Ipred=reshape(Ypred,size(I));
figure
imagesc(max(Ipred,[],3))

end

