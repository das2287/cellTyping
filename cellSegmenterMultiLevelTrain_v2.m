function [Xtrain,Ytrain,Xtest]=cellSegmenterMultiLevelTrain_v2(I,levels,numtrain)


for ds=1:length(levels)
    for z=1:size(I,3)
        I_DS{ds}(:,:,z)=imresize(imresize(I(:,:,z),levels(ds),'bilinear'),[size(I,1) size(I,2)],'bilinear');
        %         [ds z]
    end
end

Xtrain=[];
Ytrain=[];


t=0;
while t<numtrain
    randZ=randi(size(I_DS{1},3));
    imagesc(log1p(I_DS{1}(:,:,randZ)));
    title('Click on cells');
    [x,y]=getpts;
    x(end)=[];
    y(end)=[];
    coorPos=[x y repmat(randZ,[size(x,1) 1])];
    imagesc(log1p(I_DS{1}(:,:,randZ)));
    hold on
    plot(coorPos(:,1),coorPos(:,2),'r+');
    title('Click on non-cells');
    [x,y]=getpts;
    x(end)=[];
    y(end)=[];
    coorNeg=[x y repmat(randZ,[size(x,1) 1])];
    cla
    imagesc(log1p(I_DS{1}(:,:,randZ)));
    hold on
    plot(coorPos(:,1),coorPos(:,2),'r+');
    plot(coorNeg(:,1),coorNeg(:,2),'c*');
    labeledCoors=round([coorPos;coorNeg]);
    labeledCoors=labeledCoors(:,[2 1 3]);
    labels=[repmat(1,[size(coorPos,1) 1]);repmat(0,[size(coorNeg,1) 1]);];
    
    
    for i=1:length(labels)
        try
            randCoor=labeledCoors(i,:);
            for d=1:size(I_DS,2)
                tmp{d}=[I_DS{d}(randCoor(1)+1,randCoor(2),randCoor(3)) I_DS{d}(randCoor(1)+1,randCoor(2)+1,randCoor(3)) I_DS{d}(randCoor(1),randCoor(2)+1,randCoor(3)) I_DS{d}(randCoor(1)-1,randCoor(2)+1,randCoor(3)) I_DS{d}(randCoor(1)-1,randCoor(2),randCoor(3)) I_DS{d}(randCoor(1)-1,randCoor(2)-1,randCoor(3)) I_DS{d}(randCoor(1),randCoor(2)-1,randCoor(3)) I_DS{d}(randCoor(1)+1,randCoor(2)-1,randCoor(3))];
            end
            s=labels(i);
            
            
            tmp3=[];
            for d=1:size(I_DS,2)
                tmp2{d}=[];
                for g=0:7
                    tmp2{d}=[tmp2{d};circshift(tmp{d},[0 g])];
                end
                tmp3=[tmp3 tmp2{d}];
            end
            Xtrain=[Xtrain;tmp3];
            Ytrain=[Ytrain;repmat(s,[8 1])];
            t=t+1;
        catch
            warning('Cell patch out of border or some other error - going to next patch');
        end
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