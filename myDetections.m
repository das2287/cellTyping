function out=myDetections(I)


I=I(:,:,1).*imdilate(imerode(log1p(var(I,[],3))>12,strel('disk',2)),strel('disk',2));

L=bwlabel(I);
props=regionprops('table',L);

props(or(props.Area>1000,props.Area<10),:)=[];

% imagesc(I);
% hold on
% scatter(props.Centroid(:,1),props.Centroid(:,2),props.Area,'r');



out=props.Centroid;