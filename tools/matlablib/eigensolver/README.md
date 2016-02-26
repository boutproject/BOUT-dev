1. Run matlab:

    ```bash
    > module load matlab-nofonts
    > matlab
    ```

2. Run eigensolver_init.m and eigensolver_ITG.m in matlab:

    ```matlab
    >> eigensolver_init.m % assign equilibrium profiles, calculate radial basis functions, arrange mode index, evaluate matrix elements involving integrals 
    >> eigensolver_ITG.m % set up matrix and solve the eigenvalue equation for ITG
    ```

3. Run plot_eigenvalue.m and plot_eigenmode.m in matlab:

    ```matlab
    >> plot_eigenvalue.m % plot linear ITG growth rate and real frequency vs. k_theta*rho_i
    >> plot_eigenmode.m % choose toroidal mode number and plot eigenmode structure  
    ```

