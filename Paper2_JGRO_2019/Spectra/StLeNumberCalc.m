v = [0.1];
fp = [5];
d = 0.012;

Re = zeros(length(fp),length(v));
for i = 1:length(fp)
    for j = 1:length(v)
        Re(i,j) = v(j)*d/1.05E-6;     
    end
end
fprintf('Min Re #: %0.4f\n',min(min(Re)))
fprintf('Mean Re #: %0.4f\n',mean(mean(Re)))
fprintf('Max Re #: %0.4f\n',max(max(Re)))

St = zeros(length(fp),length(v));
for i = 1:length(fp)
    for j = 1:length(v)
        St(i,j) = fp(i)*d/v(j);     
    end
end
fprintf('Min St #: %0.4f\n',min(min(St)))
fprintf('Mean St #: %0.4f\n',mean(mean(St)))
fprintf('Max St #: %0.4f\n',max(max(St)))

Le = zeros(length(fp),length(v));
for i = 1:length(fp)
    for j = 1:length(v)
        Le(i,j) = v(j)/(2*pi*fp(i));     
    end
end
fprintf('Min Le: %0.4f\n',min(min(Le)))
fprintf('Mean Le: %0.4f\n',mean(mean(Le)))
fprintf('Max Le: %0.4f\n',max(max(Le)))
