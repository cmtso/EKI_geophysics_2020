function Pr=Set_prior(Grid,sigma1,sigma2,option,n_fields)
level.nu = 2.0; level.sigma = 0.5;
for idim=1:Grid.dim
    level.len(idim).mean.lim=[log(Grid.D(idim)/15),log(Grid.D(idim)/5)]; 
    if (option==0)
        level.len(idim).per.sigma= 0.0; %wheree this is zero lengthscale is not a field
    elseif (option==1)
        level.len(idim).per.sigma= 0.25; 
    end
    level.len(idim).per.nu = 2.0;
    level.len(idim).per.len{1} =log(Grid.D(idim)/7)*ones(Grid.N,1); %lengthscale of lengthscale
    level.len(idim).per.len{2} =log(Grid.D(idim)/7)*ones(Grid.N,1);   
end
Pr.sigma2=sigma2;
Pr.sigma1=sigma1;
Pr.level=level;
Pr.n_fields=n_fields;


switch Pr.n_fields
    case 3
        K(1).mean.lim=[log(1),log(1.5)]; 
        K(2).mean.lim=[log(0.05),log(0.1)]; %%%THIS IS THE BACKGROUND SIGMA (this is the value the conductivity takes outside zone 1).
        K(3).mean.lim=[log(0.01),log(0.005)]; 
        K(1).per.nu = 1.0; K(1).per.sigma =0.05;%%CHANGE HERE FOR VARIABLE SIGMA
        K(2).per.nu = 1.0; K(2).per.sigma =0.05; %%CHANGE HERE FOR VARIABLE SIGMA
        K(3).per.nu = 1.0; K(3).per.sigma =0.05; %%CHANGE HERE FOR VARIABLE SIGMA
        K(1).per.len(1).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(1).per.len(2).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(2).per.len(1).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(2).per.len(2).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(3).per.len(1).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(3).per.len(2).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 

    case 2
        K(1).mean.lim=[log(1/3000),log(1/2000)]; %%%THIS IS THE BACKGROUND SIGMA (this is the value the conductivity takes outside zone 1).
        K(2).mean.lim=[log(1/300),log(1/200)]; 
        K(1).per.nu = 1.0; K(1).per.sigma =0.05;
        K(2).per.nu = 1.0; K(2).per.sigma =0.05; %%constant for now
        
        K(1).per.len(1).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(1).per.len(2).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(2).per.len(1).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 
        K(2).per.len(2).lim=[log(Grid.D(idim)/25),log(Grid.D(idim)/5)]; 

        
end

 Pr.K=K;
