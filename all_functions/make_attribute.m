function data = make_attribute(sequence, W_IN, FHC)

% make data
filt_onesided=ones(1,W_IN);		% W_IN=21
filt_onesided(1)=1/2;
%
[data, thrash]= protein_flat_filter(sequence,[],FHC,filt_onesided);
%
for i=1:length(data(:,1))
   k2o(i,1) = entropy(data(i,1:20));
end
data(:,24)=k2o;
   
return;
