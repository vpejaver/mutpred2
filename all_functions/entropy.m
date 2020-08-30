% entropy.m

function entropy = entropy(p)

q=find(p~=0);
p=p(q);
entropy = sum(p.*log(1./p)/log(2));

return;


for i=1:size(data,1)
   e(i,1)=entropy(data(i,1:20));
end

for i=1:23
   temp=corrcoef(e',data(:,i));
   cc(i)=temp(1,2);
end
