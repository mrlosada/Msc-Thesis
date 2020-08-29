%%%%%%%%%%%%%%%%%%%%%%
% This code solves the dynamics of the atoms of a 1D array. 
% The atomic density matrix evolves under a Lindblad master equation.
% We solve the differential set of differential equations
%by using quantum jump formalism. The code is also 
%implemented to calculate the expected value of the atomic population <n_e(t)},
% the field intensity <I(t)>.Our algorithm supports a N=14 chain of 
%atoms


%Vector to store the different number of atoms
Nvec = [];

%number of trajectories
Nt = 1000;

%tMax, tVec
tMax = 1000;
vecLog = linspace(log10(0.2),log10(tMax),100);
tVec = 10.^vecLog;


%value of k and d
k = 2*pi;
d = 0.2;

%positions for the observation of I(r,t)
r = 50*d;
theta = 0;
phi = linspace(0,pi,50);



for NN = 1:length(Nvec)
    
    N = Nvec(NN);
    
    %atom positions
    z0 = (N+1)*d/2;
    z = (1:N)*d-z0;
    for i = 1:length(z)
     r_atom(:,i) = [0,0, z(i)];
    end
    
    %construction of the quantum jump operators
    gMat = zeros(N);
    for i = 1:N
        for j = 1:N
            gMat(i,j) = G(i,j,k,d);
        end
    end
    
    imG = imag(gMat);
    [Vg,gamma] = eig(imG);
    Vg = sparse(Vg);
    
    for m=1:N
        J = sparse(2^N,2^N);
        for j = 1:N
            J = Vg(j,m)*sigmage(j,N) + J;
        end
        J = sqrt(gamma(m,m))*J;
        structJ(m).matrix=J;
    end
    
    %H_eff
    Heff = sparse(2^N,2^N);
      
    for i = 1:N
        for j =1:N
            Heff = G(i,j,k,d)*sigma(i,j,N) + Heff;
        end
    end
    
  %atomic states and operators.
    kete = [1;0]; %excited state
    ketg = [0;1]; %ground state
    kete = sparse(kete);
    ketg = sparse(ketg);
    sigmae = sparse((kete')*kete); %sigmaee
    sigmag = sparse((ketg')*ketg); %sigmagg
    
    %Unitary transformation operator
    kc = pi/(2*d);
    
    for j = 1:length(kc)
        U = kete*(kete') + exp(1i*(pi+2*kc(j)*r_atom(1)))*(ketg*(ketg'));
         for i = 2:N
             U_it = kete*(kete') + exp(1i*(pi+2*kc(j)*r_atom(3,i)))*(ketg*(ketg'));
             U = kron(U,U_it);
         end
         structU_k(j).matrix = U;
    end
    
    %Intensity matrix
    for p = 1:length(phi)
        I = sparse(2^N,2^N);
        for i = 1:N
            for j = 1:N
              I = I+ conj(G_gen(r_sph(r,phi(p),theta),r_atom(:,i),k))*transpose(G_gen(r_sph(r,phi(p),theta),r_atom(:,j),k))*sigma(i,j,N);
            end
        end
        structI(p).matrix = I;
    end
    
    %excited population observable n_e = sum_n sigmae_n
    sigmaeeN = sparse(2^N,2^N);
    for i = 1:N
       sigmaeeN = sigmaeeN + sigma(i,i,N); 
    end
    
    ket0 = kete;
    %Initial full excited state |e>^{\otimes N}
    for i = 1:(N-1)
        ket0 = sparse(kron(ket0,kete));
    end
    
    %Store of propagators Ui  
    Ti = [0.001;0.01;0.1;1;10;100;1000];
    
    for i=1:length(Ti)
        structU(i,1).matrix = sparse(expm(-1i*Ti(i).*Heff));
        for p = 2:9
            disp(p)
            structU(i,p).matrix = structU(i,p-1).matrix*structU(i,1).matrix;
        end  
    end
        
    %Multidimensional array to store the expected value of the intensity
    I_ev = zeros(length(tVec),length(phi),Nt);
    I_ev_U = zeros(length(tVec),length(phi),Nt);
    I_ev_U_kc = zeros(length(kc),length(phi),Nt);

    %Dynamics
    for l = 1:Nt

        tacc = 0.2;
        ketold = ket0;

        ketJump = sparse(length(ketold),10000);
        tJump = zeros(1,1);
        tJump(1) = tacc;
        ketJump(:,1) = ketold;
        it = 2;
        
         while tacc < tVec(end)

            intm = zeros(1,m+1);

            r = rand;
            [Tau,ketold] = tau_sparse_fast(r,Ti,structU,ketold,tVec(end)-tacc);
            tacc = tacc + Tau;
            
            tJump(it) = tacc;
            if tacc > tVec(end)
                break;
            end

            mod = conj(transpose(ketold))*ketold;
            dp = 1 - mod;
            normKet = norm(ketold);
            ketold = ketold/normKet;
            
            if dp < 1e-6
                break;
            end
            
            for i = 1:m
                pm = structJ(i).matrix*ketold;
                dpm(i) = (conj(transpose(pm)))*pm;
            end

            acc = dpm;
            
            for i = 2:m
                acc(i) = acc(i)+acc(i-1);
            end

            intm = acc/acc(end);
            r2 = rand;
            j = 1;
            
            while intm(j) < r2 
                j = j +1;
            end
            
            ketold = structJ(j).matrix*ketold;
            normKet = norm(ketold);
            %normKet = conj(transpose(ketold))*ketold;
            ketold = ketold/normKet;
            ketJump(:,it) = ketold;
            it = it + 1;
         end

         tIt = 2;
         
         %Initial expected value of the intensity
         for p = 1:length(phi)

             IVal = structI(p).matrix*ketJump(:,1);
             IVal = conj(transpose(ketJump(:,1)))*IVal;
             I_ev(1,p,l) = full(IVal);
             
             if tVec(1) >= tU
                 ketAux = U*ketJump(:,1);
             else
                 ketAux = ketJump(:,1);
             end
             
             IVal_U = structI(p).matrix*ketAux;
             IVal_U = conj(transpose(ketAux))*IVal_U;
             I_ev_U(1,p,l) = full(IVal_U); 
             
         end
         
         %total excited population at the initial time
         Pe(1,l) = conj(transpose(ketJump(:,1)))*sigmaeeN*ketJump(:,1);
         Pe_U(1,l) = conj(transpose(ketAux))*sigmaeeN*ketAux;
         
         
         for i = 2:length(tJump)

             ket = ketJump(:,i-1);
             
             while tIt<=length(tVec) && tVec(tIt)<=tJump(i)
             
                 tDiff = tVec(tIt)-tVec(tIt-1);
                 
                 for pos = length(Ti):-1:1
                     q = fix(tDiff/Ti(pos));
                     r = rem(tDiff,Ti(pos));
                     if q ~= 0
                        ket = structU(pos,q).matrix*ket;
                     end
                     tDiff = r;
                 end

                 ket_U = structU_k(1).matrix*ket;

                 %Expected value of the intensity for the transformed
                 %state.
                 for p = 1:length(phi)
                   I_ev_U(tIt,p,l) = conj(transpose(ket_U))*structI(p).matrix*ket_U;
                 end
                 
                 %Expected value of the intensity for the original state
                 for p = 1:length(phi)
                   I_ev(tIt,p,l) = conj(transpose(ket))*structI(p).matrix*ket;
                 end
                 
                 %total excited population
                 Pe(tIt,l) = conj(transpose(ket))*sigmaeeN*ket;
                 Pe_U(tIt,l) = conj(transpose(ket_U))*sigmaeeN*ket_U;
                 tIt = tIt+1;
             end
             
         end
         
    end
   teP(:,NN) = mean(Pe,2);
   teP_U(:,NN) = mean(Pe_U,2);
   ieV(:,:,NN) = mean(I_ev,3);
   ieV_U(:,:,NN) = mean(I_ev_U,3);
   ieV_U_kc(:,:,NN) = mean(I_ev_U_kc,3);
   
end


for NN = 1:length(Nvec)
    
    for i = 1:length(tVec)
        ieVNorm(i,:,NN) = ieV(i,:,NN)/norm(ieV(i,:,NN));
        tNorm(i) = norm(ieV(i,:,NN));
        ieV_U_Norm(i,:,NN) = ieV_U(i,:,NN)/norm(ieV_U(i,:,NN));
    end
end

%Integration over the solid angle of the angular distribution of the 
%intensity to obtain the total intensity I_total as a function of time.

for NN = 1:length(Nvec)
    I_total = zeros(length(tVec),length(Nvec));
    for i = 1:length(tVec)
        for j = 1:length(phi)
            I_total(i,NN) = ieV_U(i,j,NN).*sin(phi(j)) + I_total(i,NN);
        end
        I_total(i,NN) = I_total(i,NN)*2*pi/length(phi);
    end
    
end

save('teP.mat','teP')
save('ieV.mat','ieV')
save('ieV_U.mat','ieV_U')
save('ieVNorm.mat','ieVNorm')
save('ieV_U_Norm.mat','ieV_U_Norm')
save('Int.mat','Int')
