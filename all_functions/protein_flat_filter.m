% protein_flat_filter.m     02/25/00

% using results from protein_read makes flat file of protein data
% Format:   (20 filtered attributes of AAs;
%            3 filtered attributes for flexibility, hydropathy, coordinate;
%            target (order/disorder);
%            number of protein; position of AA in a given protein)


function [order, disorder]= protein_flat_filter(sequence,disorder_pos,FHC,filt_onesided)

l=length(sequence);
order=[];
disorder=[];
if isempty(disorder_pos)
   disorder_pos{l}=[];
end
IO=1;
ID=1;
for i=1:l
   matrix_protein=[];
   temp=double(sequence{i});
   u=unique(temp);
   for j=1:length(u)
     f=find(FHC(1,:)==u(j));
      q=find(temp==u(j));
      matrix_protein(q,f)=1;
      matrix_protein(q,21:23)=repmat(FHC(2:4,f)',length(q),1);     
    end
    if isempty(disorder_pos{i})==1
      matrix_protein(:,24)=zeros(length(sequence{i}),1);
      % positional data, just for a reference...
      matrix_protein(:,25)=IO*ones(length(sequence{i}),1);
      IO=IO+1;
      matrix_protein(:,26)=(1:length(sequence{i}))';
      % filtering...
      matrix_protein(:,1:23)=filter_twosided(matrix_protein(:,1:23),filt_onesided);
      % add to data...
      order=[order; matrix_protein];     
    else
      l_dis=length(disorder_pos{i});
      dis_pos=[];
      for k=1:l_dis/2
        dis_pos=cat(1,dis_pos,(disorder_pos{i}(2*k-1):disorder_pos{i}(2*k))');
      end
      matrix_protein(dis_pos,24)=1;
      % positional data, just for a reference...
      matrix_protein(:,25)=ID*ones(length(sequence{i}),1);
      ID=ID+1;
      matrix_protein(:,26)=(1:length(sequence{i}))';
      % filtering...
      matrix_protein(:,1:23)=filter_twosided(matrix_protein(:,1:23),filt_onesided);
      % add to data...
      disorder=[disorder; matrix_protein];      
    end    
end

return;

% to have mapping from 20 AAs to 3 important attributes
% Flexibility, Hydropathy, Coordinate
FHC(1,:)=double('ACDEFGHIKLMNPQRSTVWY')

% for F
tmp(1,:)=double('WCFIYVLHMATRGQSNPDEK')
tmp(2,:)=[.94 .906 .915 .927 .929 .931 .935 .95 .952 .984 .997 1.008 1.031 1.037 1.046 1.048 1.049 1.068 1.094 1.102];
[a,b]=sort(tmp(1,:)',1);
FHC(2,:)=tmp(2,b);

% for H
tmp(1,:)=double('IVLFCMAGTWSYPHEQDNKR')
tmp(2,:)=[4.5 4.2 3.8 2.8 2.5 1.9 1.8 -0.4 -0.7 -0.9 -0.8 -1.3 -1.6 -3.2 -3.5 -3.5 -3.5 -3.5 -3.9 -4.5]
[a,b]=sort(tmp(1,:)',1);
FHC(3,:)=tmp(2,b);

% for C
tmp(1,:)=double('WCFIYVLHMATRGQSNPDEK')
tmp(2,:)=[5.81 7.49 5.99 5.92 6.13 5.60 6.00 6.11 4.89 4.97 5.00 4.75 4.72 4.82 4.31 5.02 3.71 3.99 3.97 4.47]
[a,b]=sort(tmp(1,:)',1);
FHC(4,:)=tmp(2,b)

save FHC FHC


%%%%%%%%%%%%%%%%%%%%%%%%%
% choose one-sided filter

F=10;
filt_onesided=2*ones(F,1)/(2*F-1);
filt_onesided(1)=filt_onesided(1)/2;
filt_onesided=[0.5 .9 .8 .7 .6 .5 .4 .3 .2 .1];
data=a;
