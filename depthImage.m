function I=depthImage(V)

colors=hsv(size(V,3));
I=zeros(size(V,1),size(V,2),3);
for ch=1:3
    I(:,:,ch)=sum(V.*reshape(colors(:,ch),[1 1 size(V,3)]),3);
end