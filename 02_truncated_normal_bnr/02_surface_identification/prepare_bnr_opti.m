clear all
clc

load('optimal_brn_25reps.mat');

[nn,rr] = meshgrid(double(n),r_values);

n = double(n);
r = r_values;
rbn_MC_data = squeeze(mean(brn_matrix,1));

figure(1)
surf(nn,rr,rbn_MC_data);
xlabel('\it{n}');
ylabel('\it{r}');
zlabel('\it{b(n,r)}');

diary bnr_data.txt

[i1,i2] = size(nn);

fprintf('nn = [ ...\n');
for i =1:i1
  for j = 1:i2
    if j<i2
       fprintf('%i, ',nn(i,j));
    else
    fprintf('%i',nn(i,j));
    end

   if j./10==round(j./10)
     fprintf('... \n');
   end

  end
  fprintf('; ... \n');
end
%
fprintf(']; \n\n');

[i1,i2] = size(rr);

fprintf('rr = [ ...\n');
for i =1:i1
  for j = 1:i2
    if j<i2
       fprintf('%1.2f, ',rr(i,j));
    else
    fprintf('%1.2f',rr(i,j));
    end

   if j./10==round(j./10)
     fprintf('... \n');
   end

  end
  fprintf('; ... \n');
end
%
fprintf(']; \n\n');

[i1,i2] = size(rbn_MC_data);
fprintf('bnr = [ ...\n');
for i =1:i1
  for j = 1:i2
    if j<i2
       fprintf('%1.5f, ',rbn_MC_data(i,j));
    else
    fprintf('%1.5f',rbn_MC_data(i,j));
    end

   if j./10==round(j./10)
     fprintf('... \n');
   end

  end
  fprintf('; ... \n');
end
%
fprintf(']; \n\n');


bnr_std = squeeze(std(brn_matrix,1));

figure(2)
surf(nn,rr,bnr_std);
xlabel('\it{n}');
ylabel('\it{r}');
zlabel('\it{s(b(n,r))}');


fprintf('bnr_std = [ ...\n');
for i =1:i1
  for j = 1:i2
    if j<i2
       fprintf('%1.5f, ',bnr_std(i,j));
    else
    fprintf('%1.5f',bnr_std(i,j));
    end

   if j./10==round(j./10)
     fprintf('... \n');
   end

  end
  fprintf('; ... \n');
end
%
fprintf(']; \n\n');

diary off

