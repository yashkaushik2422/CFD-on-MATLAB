clc
clear 
close all


% Initialize parameters
imAll = 100; % Number of cells in x-direction
jmAll = 100; % Number of cells in y-direction
ng = 2;      % Number of ghost cells (both sides, even)
delta = 1;   % Grid spacing
rho = 10;    % Density
CFL = 0.1;   % Courant–Friedrichs–Lewy number.
D = 5*(10.^-1); %Diffusivity

% Calculate indices and grid
    imTotal = imAll + ng;
    ifi = ng / 2 + 1;
    ila = ifi + imAll - 1;
    ifim = ifi - 1;
    ilam = ila - 1;
    ifip = ifi + 1;
    ilap = ila + 1;
    
    jmTotal = jmAll + ng;
    jfi = ng / 2 + 1;
    jla = jfi + jmAll - 1;
    jfim = jfi - 1;
    jlam = jla - 1;
    jfip = jfi + 1;
    jlap = jla + 1;
    
    dx = delta;
    dy = delta;
    
    x = linspace(dx / 2, (imAll - 1 / 2) * dx, imAll);
    y = linspace(dy / 2, (jmAll - 1 / 2) * dy, jmAll);



    [X, Y] = meshgrid(x, y);

  % Initialize velocities and scalar

    phiU = zeros(imTotal, jmTotal);
    phiV = zeros(imTotal, jmTotal);
    phiU_New = zeros(imTotal, jmTotal);
    phiV_New = zeros(imTotal, jmTotal);
    
    %U and V with initial values
    U = sin(2*pi ./ imAll * X) .* cos(2*pi .* Y./ jmAll); %Taylor-Green-Vortex
    V = cos(2*pi ./ jmAll * X) .* sin(2*pi .* Y./ jmAll);

    %test velocities
    %U = ones(imAll,imAll);
    %V = zeros(imAll,imAll);
    
    % Compute phiU and phiV (momentum)
    phiU(ifi:ila, jfi:jla) = rho .* U;           %sin(pi*Y/100);
    phiV(ifi:ila, jfi:jla) = rho .* V;

    phiU(ifim,:) = phiU(ila,:); 
    phiU(ilap,:) = phiU(ifi,:);
    phiU(:,jfim) = phiU(:,jla); 
    phiU(:,jlap) = phiU(:,jfi);

    
    %initialize time-step
    dt = CFL .* delta./max(max(abs(U)));

    % For n-iterations in time-stepping

for n = 1 : 100000

    %For one iteration
    FX = phiU;
    FY = phiV;
    
    %Convective flux by CDS
    fluxConX =  ( FX(ifip:ilap,jfi:jla) .* U - FX(ifim:ilam,jfi:jla) .* U );  %*0.5./dx *
    fluxConY =  ( FY(ifi:ila,jfip:jlap) .* V - FY(ifi:ila,jfim:jlam) .* V );  % 0.5./dy *
    fluxConNet = (fluxConX + fluxConY);

    % Diffusive flux by CDS
    fluxDifX = ( (FX(ifim:ilam,jfi:jla) - 2.*FX(ifi:ila,jfi:jla) + FX(ifip:ilap,jfi:jla)) ); % D./(dy).^2 
    fluxDifY = ( (FY(ifi:ila,jfip:jlap) - 2.*FY(ifi:ila,jfi:jla) + FY(ifi:ila,jfim:jlam)) );
    fluxDifNet = (fluxDifY + fluxDifY);

    
    % Update phi as per convection diffusion

    phiU_New(ifi:ila,jfi:jla) = phiU(ifi:ila,jfi:jla) + dt.*( D./(dx).^2.*(fluxDifNet) - 0.5./dx .*(fluxConNet) );%multiplation with scalar here instead of the function file
    
    phiV_New(ifi:ila,jfi:jla) = phiV(ifi:ila,jfi:jla) + dt.*( D./(dy).^2.*(fluxDifNet) - 0.5./dy .*(fluxConNet) ); %multiplation with scalar here instead of the function file

    
    %divergence prediction of scalar
    divPredU(ifi:ila,jfi:jla) =  (0.5./dx)*(phiU(ifip:ilap,jfi:jla) - phiU(ifim:ilam,jfi:jla));
    divPredV(ifi:ila,jfi:jla) =  (0.5./dy)*(phiV(ifi:ila,jfip:jlap) - phiV(ifi:ila,jfim:jlam));
    divPredMom =  divPredV + divPredU; 

    %Poisson solver   (P''_x = [RhoU'_x]_t) for U direction

        P = zeros(imTotal, imTotal);
        residual = inf;  
        tolerance = 1e-6; % Convergence tolerance
        max_iterations = 1000; % Maximum number of iterations to prevent infinite loop
        iteration = 0; % Initialize iteration counter
        
        while residual > tolerance && iteration < max_iterations
            iteration = iteration + 1;
            P_old = P; 
        
            P(ifi:ila,jfi:jla) = 0.25 * (P(ifim:ilam,jfi:jla) + P(ifip:ilap,jfi:jla) + P(ifi:ila,jfim:jlam) + P(ifi:ila,jfip:jlap)  - divPredMom(ifi:ila,jfi:jla) * (dx^2) ./ dt);
        
            residual = norm(P - P_old, 'fro'); % Frobenius norm for the residual
        end

    %Correction of Momentum
    phiU(ifi:ila,jfi:jla) = phiU(ifi:ila,jfi:jla) - dt.*( P(ifip:ilap,jfi:jla) - P(ifim:ilam,jfi:jla) )./(2*dx); 
    phiV(ifi:ila,jfi:jla) = phiV(ifi:ila,jfi:jla) - dt.*( P(ifi:ila,jfip:jlap) - P(ifi:ila,jfim:jlam) )./(2*dy); 
    
    %Velocity Correction [Ask if necessary]
    U = phiU(ifi:ila,jfi:jla)./rho; 
    V = phiV(ifi:ila,jfi:jla)./rho;

    %Velocity correction with Center2Surface & Mom2Vel
    phiU_surface = (phiU(ifim:ilam,ifim:ilam)+ phiU(ifip:ilap,ifip:ilap))./2;
    phiV_surface = (phiV(ifim:ilam,ifim:ilam)+ phiV(ifip:ilap,ifip:ilap))./2; 
    
    U_surface = phiU_surface ./rho;
    V_surface = phiV_surface ./rho;


    % Apply periodic boundaries (ghost-cells)
    phiU_New(ifim,:) = phiU_New(ila,:); 
    phiU_New(ilap,:) = phiU_New(ifi,:);
    phiU_New(:,jfim) = phiU_New(:,jla); 
    phiU_New(:,jlap) = phiU_New(:,jfi);

    phiV_New(ifim,:) = phiV_New(ila,:); 
    phiV_New(ilap,:) = phiV_New(ifi,:);
    phiV_New(:,jfim) = phiV_New(:,jla); 
    phiV_New(:,jlap) = phiV_New(:,jfi);
    
    %Move one time step ahead by updating Phi
    phiU = phiU_New;
    phiV = phiV_New;


     %visualization
    % if (mod(n,10) == 0)
    %    imagesc(phiU);
    % end
    % pause(0.02);

    % visualization
    if (mod(n,10) == 0)
        imagesc(phiV);
    end
    pause(0.02);


end















