function [I_Z,I_S,I_1st_moment,I_2nd_moment]=spatial_zscore(I,s)
tic;I_1st_moment=imgaussfilt3(I,s);
I_2nd_moment=imgaussfilt3(I.^2,s);
I_S=sqrt(max(I_2nd_moment-I_1st_moment.^2,eps));
I_Z=(I-I_1st_moment)./I_S;toc
end