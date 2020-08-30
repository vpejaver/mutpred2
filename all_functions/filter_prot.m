% take needed statistics regarding predictions

function pfilt = filter_prot(pred,prot_no,filt_onesided)

%
% Calculating what is the number of proteins in the database 
% with chance of long disorder
%
pfilt=[];
dd=diff(prot_no);
dd(1)=1;
q=find(dd~=0);
q(1)=0;
q(length(q)+1)=length(dd)+1;
III=1;
for i=1:length(q)-1
   if mod(i,50)==0
      50*III
      psave{III}=pfilt;
      pfilt=[];
      III=III+1;
   end   
   % predictions for protein i
   p=pred(q(i)+1:q(i+1));
   clear pf
   pf(:,1)=filter_error(p,filt_onesided,0);
   pfilt=[pfilt;pf];
end
psave{III}=pfilt;

if III>1
   pfilt=[];
   L=length(psave);
   for i=1:L
      pfilt=[pfilt; psave{i}];
   end
end

return  


a=[];
for i=1:length(dis_all)
   a=[a;dis_all{i}(:,25)];
end
