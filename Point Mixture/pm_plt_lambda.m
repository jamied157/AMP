function [] = pm_plt_lambda(e,lambda)
n = 2000;
y = zeros(n,length(lambda));
gamma = zeros(1,length(lambda));
x = zeros(n,length(lambda));
z = zeros(1,length(lambda));

for i = 1:length(lambda)
    [A,x(:,i)] = ex_pmatrix(n,lambda(i),e);
    [y(:,i), gamma(i)] = pm_factor(A,e,lambda(i));
    z(i) = abs(dot(y(:,i),x(:,i))/(norm(y(:,i))*norm(x(:,i))));
end
scatter(lambda,z.^2,'ko','MarkerFaceColor','r')
hold on
plot(lambda,gamma./lambda.^2,'b','LineWidth',2)
legend('empirical overlap','theoretical overlap','Location','east')
title(['\epsilon = ' ,num2str(e)])
xlabel('\lambda')
ylabel('Overlap')
ylim([0,1])
hold off
end