q1=load('Q1'); q2=load('Q2'); q3=load('Q3'); q4=load('Q4'); q5=load('Q5');
[t,n]=size(q1);

txy=zeros(t,n);
tyz=zeros(t,n);
S=zeros(t,n);
for i=1:t
    for j=1:n
        a1=q1(i,j); a2=q2(i,j); a3=q3(i,j); a4=q4(i,j); a5=q5(i,j);
        Q=[a1,a2,a3; a2,a4,a5; a3,a5,-a1-a4];
        [v lam]=eig(Q);
        lam=diag(lam);
        [s k]=max(lam);
        v0=v(:,k);
        txy(i,j)=atan2(v0(2),v0(1))/pi*180;
        if txy(i,j) < 0
            txy(i,j) = txy(i,j) + 180;
        end
        tyz(i,j)=atan2(v0(3),v0(2))/pi*180;
        if tyz(i,j) < 0
            tyz(i,j) = tyz(i,j) + 180;
        end
        
        S(i,j)=s;
    end
end

save('txy','-ascii','txy');
save('tyz','-ascii','tyz');
save('S','-ascii','S');
