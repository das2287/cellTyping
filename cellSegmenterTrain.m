function [Xtrain,Ytrain,Xtest]=cellSegmenterTrain(I,numtrain)
patchRadius=30;
[x,y,z]=meshgrid(1:size(I,1),1:size(I,2),1:size(I,3));

x=permute(x,[2 1 3]);
y=permute(y,[2 1 3]);
z=permute(z,[2 1 3]);
coor0=[x(:) y(:) z(:)];


Xtrain=[];
Ytrain=[];


for q=0:(max(log1p(I(:)))-1)
    t=0;
    coor=coor0(log1p(I)>q,:);
    while t<numtrain
        clear patch
        try
            randCoor=coor(randi(size(coor,1)),:); %% add more criteria here if needed
            tmp=[I(randCoor(1)+1,randCoor(2),randCoor(3)) I(randCoor(1)+1,randCoor(2)+1,randCoor(3)) I(randCoor(1),randCoor(2)+1,randCoor(3)) I(randCoor(1)-1,randCoor(2)+1,randCoor(3)) I(randCoor(1)-1,randCoor(2),randCoor(3)) I(randCoor(1)-1,randCoor(2)-1,randCoor(3)) I(randCoor(1),randCoor(2)-1,randCoor(3)) I(randCoor(1)+1,randCoor(2)-1,randCoor(3))];
            patch=I(randCoor(1)-patchRadius:randCoor(1)+patchRadius,randCoor(2)-patchRadius:randCoor(2)+patchRadius,randCoor(3));
            imagesc(patch);title('Is this a cell?');
            hold on
            plot(patchRadius,patchRadius,'r+','MarkerSize',20);
            s=input('Is this a cell?(1/0)');
            if isempty(s)
                s=0;
            end
            if s~=1
                s=0;
            end
            t=t+1;
            tt=0;
            tmp2=[];
            for g=0:7
                tmp2=[tmp2;circshift(tmp,[0 g])];
            end
            Xtrain=[Xtrain;tmp2];
            Ytrain=[Ytrain;repmat(s,[8 1])];
        catch
            warning('Cell patch out of border or some other error - going to next patch');
        end
    end
    
end
shifts=[1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1];
Xtest=zeros(numel(I),size(shifts,1));
for s=1:size(shifts,1)
    Ishift=imtranslate(I,shifts(s,:));
    Xtest(:,s)=Ishift(:);
end



end