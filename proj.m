% projection operator. JJ Z
function xp=proj(a,k)
if max(abs(a))==0
    xp=0;
else
    [~,index1]=sort(abs(a),'descend');
    tka=zeros(size(a));
    for j=1:k
        tka(index1(j))=a(index1(j));
    end
    xp=tka/norm(tka);
end
end