clear all

load('Results_std_brn.mat');

[nn,rr] = meshgrid(nT,rT);

figure(1)
surf(nn,rr,sm_arT);
xlabel('\it{n}');
ylabel('\it{r}');
zlabel('\it{s-a(r)}');

figure(2)
surf(nn,rr,sm_ar_bnrT);
xlabel('\it{n}');
ylabel('\it{r}');
zlabel('\it{s-(b(n,r),a(r))}');


