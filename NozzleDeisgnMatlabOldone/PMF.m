function [ M,nu,mu ] = PMF(G,M,nu,mu)

Gp=G+1;
Gm=G-1;

% for known M
if M~=0;

    nu = sqrt(Gp/Gm).*atand(sqrt(Gm*(M.^2-1)/Gp))-atand(sqrt(M.^2-1));

    mu = asind(1./M);
    

% for known nu
elseif norm(nu)~=0;
    
    % Find M
        
    %Nu = @(Mg)sqrt(Gp/Gm)*atand(sqrt(Gm*(Mg.^2-1)/Gp))-atand(sqrt(Mg.^2-1))-nu;
    
    for i=1:length(nu(1,:))
        for j = 1:length(nu(:,1))
            M(j,i) = fzero(@(Mg)sqrt(Gp/Gm)*atand(sqrt(Gm*(Mg.^2-1)/Gp))...
                -atand(sqrt(Mg.^2-1))-nu(j,i),[1 100]);
        end
    end

    mu = asind(1./M);
    
    
% for known mu
elseif mu~=0;
    
    M=1./sind(mu);
    
    nu=sqrt(Gp/Gm)*atand(sqrt(Gm*(M.^2-1)/Gp))-atand(sqrt(M.^2-1));
    
end