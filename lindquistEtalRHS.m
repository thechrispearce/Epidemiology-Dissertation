function dy=lindquistEtalRHS(y,beta,gamma,M)

% calculate the right hand side of the Lindquist etal ODE model for the
% SIS-CM epidemic.

tfm=@(x,i,j)lindquistEtalStateTfm(M,x,i,j);

%first calculate the values of G and H, the big fractions

Gnum=0;
Gden=0;
Hnum=0;
Hden=0;

for j=0:M
    for k=0:M-j
        Gnum=Gnum+j*k*y(tfm(0,j,k));
        Gden=Gden+j*y(tfm(0,j,k));
        Hnum=Hnum+k^2*y(tfm(0,j,k));
        Hden=Hden+j*y(tfm(1,j,k));
    end
end

G=beta*Gnum/Gden;
H=beta*Hnum/Hden;

dy=zeros((M+1)*(M+2),1);

for s=0:M
    for i=0:M-s
        %first do S'_{si}
        dy(tfm(0,s,i)) = -y(tfm(0,s,i))*((beta+gamma)*i+G*s) ...
                         +gamma*y(tfm(1,s,i));
        if i<M && s>0
            dy(tfm(0,s,i)) = dy(tfm(0,s,i)) + gamma*(i+1)*y(tfm(0,s-1,i+1));
        end
        
        if i>0 && s<M
            dy(tfm(0,s,i))= dy(tfm(0,s,i)) + G*(s+1)*y(tfm(0,s+1,i-1));
        end
        
        %now do I'_{si}
        dy(tfm(1,s,i)) = beta*i*y(tfm(0,s,i)) ...
                        - y(tfm(1,s,i))*(gamma*(i+1)+H*s);
        if s>0 && i<M
            dy(tfm(1,s,i)) = dy(tfm(1,s,i)) + gamma*(i+1)*y(tfm(1,s-1,i+1));
        end
        
        if s<M && i>0
            dy(tfm(1,s,i)) = dy(tfm(1,s,i)) + H*(s+1)*y(tfm(1,s+1,i-1));
        end
    end
end