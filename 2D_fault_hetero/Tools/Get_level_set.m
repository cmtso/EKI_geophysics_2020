function U=Get_level_set(Grid,Pr,L_means,L_per,RN)

    pri.len{1}=L_per.x+L_means(1);
    pri.len{2}=L_per.y+L_means(2);
    pri.sigma=Pr.level.sigma; pri.nu=Pr.level.nu;
    U=reshape(grf2D(Grid, pri,RN),Grid.n(1),Grid.n(2));
