clear all
clc
p = 2.34;

k = 5e6;      % number of tries
nT = 2:20;     % check points

sqr = @(x) x.*x;

targetF = @(x) sum(sqr(x-p));

std1  = @(x) sqrt(sum(sqr(x-mean(x)))./(size(x,1)-1));      % n-1
std15 = @(x) sqrt(sum(sqr(x-mean(x)))./(size(x,1)-1.5));    % n-1.5
stdV  = @(x) sqrt(sum(sqr(x-mean(x)))./(size(x,1)-1.5+1./(8.*(size(x,1)-1.5)))); % Wikipedia
stdo1  = @(x) sqrt(sum(sqr(x-mean(x)))./(size(x,1)-(1.521.*size(x,1)-1.6415)./(1.014.*size(x,1)-1)) ); % optim1


stdT   = [];
std1T  = [];
std15T = [];
stdVT  = [];
stdo1T = [];
stdo2T = [];

fprintf('Modelling started for k= %1.2e repetitions:\n',k);
for ni = 1:numel(nT)
  fprintf(' ');
end
%
printf('*\n');

start_time = time();

for ni = 1:numel(nT)
  n = nT(ni);

  T = randn(n,k).*p;

  stdT   = [stdT;   mean(std(T))  ];  % implemented std
  std1T  = [std1T;  mean(std1(T)) ];  % std1 - verification
  std15T = [std15T; mean(std15(T))];  % std15 - verification
  stdVT  = [stdVT;  mean(stdV(T)) ];  % stdV - verification
  stdo1T = [stdo1T; mean(stdo1(T))];  % stdo1 - verification
  stdo2T = [stdo2T; mean(std_opti(T))];  % std_opti - verification

  fprintf('.');
end
%

elapsed_time = time() - start_time;
printf('\nModelling done in %.2f seconds\n\n', elapsed_time);

printf('Parameters quality:\n');
printf('original:                        %1.5e\n', targetF(stdT)  );
printf('n-1:                             %1.5e\n', targetF(std1T) );
printf('n-1.5:                           %1.5e\n', targetF(std15T));
printf('n-1.5+1/(8*(n-1)):               %1.5e\n', targetF(stdVT) );
printf('n-(1.521*n-1.6415)/(1.014*n-1)): %1.5e\n', targetF(stdo1T) );
printf('n-b(n):                          %1.5e\n', targetF(stdo2T) );
printf('\n\n');


figure (1)
plot (nT, stdT,  '-*r', ...
      nT, std1T, '-or', ...
      nT, std15T,'-om', ...
      nT, stdVT, '-ob', ...
      nT, stdo1T,'-sb', ...
      nT, stdo2T,'-sr');
xlim([2,max(nT)]);
% ylim([0 1.1]);
grid;
xlabel('n');
ylabel('std dev estimator');
legend({'stdT','std1T','std15T','stdVT','stdo1T','std opti'});

figure (2)
semilogy(nT, abs(stdT-p),  '-*r', ...
         nT, abs(std1T-p), '-or', ...
         nT, abs(std15T-p),'-om', ...
         nT, abs(stdVT-p), '-ob', ...
         nT, abs(stdo1T-p),'-sb', ...
         nT, abs(stdo2T-p),'-sr');
xlim([2,max(nT)]);
% ylim([0 1.1]);
grid;
xlabel('n');
ylabel('std dev estimator error');
legend({'stdT','std1T','std15T','stdVT','stdo1T','std opti'});



