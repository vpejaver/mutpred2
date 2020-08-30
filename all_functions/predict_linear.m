% predict_linear.m
%
% INPUT:
%	models
%  sequence
% OUTPUT:
%  pfilt - predicted sequence

function pfilt = predict_linear(sequence, models)

global CURRDIR

% make data
load(strcat(CURRDIR, filesep, 'all_models', filesep, 'FHC.mat'));
W_IN = 21;
filt_onesided=ones(1,W_IN);		% W_IN=21
filt_onesided(1)=1/2;
%
s{1}=sequence;
[data, thrash]= protein_flat_filter(s,[],FHC,filt_onesided);
%
for i=1:length(data(:,1))
   k2o(i,1) = entropy(data(i,1:20));
end
data(:,22)=k2o;

W_OUT = 21;

filt_onesided = ones(1,W_OUT);		
filt_onesided(1) = 1/2;

lmod=length(models);
Attributes = [1:3 5 7:20 21 22];
for j=1:lmod
   p(:,j)=[ones(length(data(:,1)),1) data(:,Attributes)]*models{j}.beta;      
	pfilt(:,j)=filter_prot(p(:,j), ones(length(p(:,j)),1), filt_onesided);
end
   
return;

