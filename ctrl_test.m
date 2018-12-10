x0=[1;0;1;0];
f = @(t,x)A*x+B*K*(x-xd);
[t,y]=ode45(f, [0 5], x0);
plot(y(:,1),y(:,2))