clear all

if pkg_ready('statistics')
  fprintf('Pkg statistics loaded \n');
else
  fprintf('ERROR:" Pkg statistics not installed \n');
  return;
end

k = 1e6;

nT = 2:20;
rT = [0.3 0.5 1 2 5];

s1 = zeros(k,1);
s2 = s1;

sm_ar_bnr = nT.*0;
sm_ar = nT.*0;

sm_ar_bnrT = [];
sm_arT = [];

tic;

for ri_ = 1:numel(rT)
  r = rT(ri_);
  for ni_ = 1:numel(nT)
    n = nT(ni_);

    for i=1:k

      x = zeros(n,1);
      for j = 1:n
        x(j) = trandn(-1.*r,r);
      end

      s1(i) = std_bnr_opti(x,r);
      s2(i) = std_bnr(x,r);

    end
    %
    fprintf('Position: n = %i, r = %1.3f, k = %i\n',n,r,k);
    fprintf('s_ar_bnr = %1.6f \n',mean(s1));
    fprintf('s_ar     = %1.6f \n',mean(s2));
    fprintf('\n');

    sm_ar_bnr(ni_) = mean(s1);
    sm_ar(ni_) = mean(s2);

  end
  %
  sm_ar_bnrT = [sm_ar_bnrT; sm_ar_bnr];
  sm_arT = [sm_arT; sm_ar];

end
%

toc

save -v7 Results_std_brn.mat sm_ar_bnrT sm_arT k nT rT





