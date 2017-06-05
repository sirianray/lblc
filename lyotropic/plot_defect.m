dpos1=load('./trial1n/dpos');
dpos2=load('./trial2/dpos');
dpos3=load('./trial3/dpos');

figure(2);clf;

plot(dpos3(1:2:end,1),'bo-','markerfacecolor','b');hold on
plot(dpos3(1:2:end,4),'bo-');hold on
plot(dpos1(1:2:end,1),'rs-','markerfacecolor','r');hold on
plot(dpos1(1:2:end,4),'rs-');hold on
plot(dpos2(1:2:end,1),'gd-','markerfacecolor','g');hold on
plot(dpos2(1:2:end,4),'gd-');hold on
line([0 200],[dpos3(end,1) dpos3(end,1)],'color','b','linestyle','--')
line([0 200],[dpos1(end,1) dpos1(end,1)],'color','r','linestyle','--')
line([0 200],[dpos2(end,1) dpos2(end,1)],'color','g','linestyle','--')
hold off

ldg=legend('+1/2, K_{11}/K_{33}=0.3','-1/2, K_{11}/K_{33}=0.3',...
    '+1/2, K_{11}/K_{33}=1.0','-1/2, K_{11}/K_{33}=1.0',...
    '+1/2, K_{11}/K_{33}=3.3','-1/2, K_{11}/K_{33}=3.3','location','southeast');
set(ldg,'box','off')
xlabel('simulation time','fontsize',20)
ylabel('x (\xi)','fontsize',20)