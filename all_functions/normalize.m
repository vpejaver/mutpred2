% normalize.m , 03/20/98

% if meanvv, stdvv are empty, normalizes data by columns
% and returns mean and std values
% else uses meanvv, stdvv to rescale data

function [meanv,stdv,datanew]=normalize(data,meanvv,stdvv)

a=size(data);
e=ones(a(1),1);
if length(meanvv)==0
   meanv=mean(data);
   stdv=std(data);
   % nabudzeno za ovaj problem
   q=find(stdv<0.000000001);
   stdv(q)=1000;
else
   meanv=meanvv;
   stdv=stdvv;
end

for i=1:a(2)
   x=data(:,i)-e*meanv(i);
   datanew(:,i)=x/stdv(i);
end

return;