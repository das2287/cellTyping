function [Xtrain,Ytrain,Xtest]=cellSegmenterMultiLevelTrain(I,levels,numtrain)


for ds=1:length(levels)
    for z=1:size(I,3)
        I_DS{ds}(:,:,z)=imresize(imresize(I(:,:,z),levels(ds),'bilinear'),[size(I,1) size(I,2)],'bilinear');
%         [ds z]
    end
end


patchRadius=30;
[x,y,z]=meshgrid(1:size(I_DS{1},1),1:size(I_DS{1},2),1:size(I_DS{1},3));

x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);
coor0=[x(:) y(:) z(:)];


Xtrain=[];
Ytrain=[];


for q=0:0.5:(max(log1p(I(:)))-1)
    t=0;
    coor=coor0(log1p(I)>q,:);
    while t<numtrain
        clear patch tmp
        try
            randCoor=coor(randi(size(coor,1)),:); %% add more criteria here if needed
            for d=1:size(I_DS,2)
                tmp{d}=[I_DS{d}(randCoor(1)+1,randCoor(2),randCoor(3)) I_DS{d}(randCoor(1)+1,randCoor(2)+1,randCoor(3)) I_DS{d}(randCoor(1),randCoor(2)+1,randCoor(3)) I_DS{d}(randCoor(1)-1,randCoor(2)+1,randCoor(3)) I_DS{d}(randCoor(1)-1,randCoor(2),randCoor(3)) I_DS{d}(randCoor(1)-1,randCoor(2)-1,randCoor(3)) I_DS{d}(randCoor(1),randCoor(2)-1,randCoor(3)) I_DS{d}(randCoor(1)+1,randCoor(2)-1,randCoor(3))];
                patch{d}=I_DS{d}(randCoor(1)-patchRadius:randCoor(1)+patchRadius,randCoor(2)-patchRadius:randCoor(2)+patchRadius,randCoor(3));
                subplot(1,size(I_DS,2),d)
                imagesc(patch{d});title('Is this a cell?');
                hold on
                plot(patchRadius,patchRadius,'r+','MarkerSize',20);
            end
            s=input('Is this a cell?(1/0)');
            if isempty(s)
                s=0;
            end
            if s~=1
                s=0;
            end
            
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