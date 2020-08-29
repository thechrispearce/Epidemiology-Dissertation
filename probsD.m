function varargout=probsD(type,a,i,m)

% y=probsD(type,a,i,m) returns the probability of the degree distribution
% taking the value i. type is a string saying what kind of rv D is and a is
% a vector of parameters for it to use. m modifies the distribution :
% m=0 D as-is,
% m=1 \tilde{D}-1
% m=2 \tilde{D},
%
% current choices are poi(lambda), const(alpha), geom(p), zeta(r),
% given(vector of P(D=i) i=0,1,2,...), heavy([k,a]), heavyC([kappa,a]),
% shifted({type,param,shift})

%global stirling1;

y=0;

if (i<0) || (m==2 && i==0)
    varargout{1}=y;
    
    %if the call is expecting paramD to be modified then say 'no'
    if nargout>1
        error('sirl:probsD:noGood','You''re expecting paramD to be modified but the call to probsD trivially returns 0 so doesn''t get that far')
    end
    return
end


switch type
    case 'heavyC' % heavy-tailed with exponential cutoff
        if m==0
            if length(a)>2 && a(3)>0 %if constant already calculated
                c=a(3);              %  use it
            else                     %else calculate it
                if abs(a(2)-round(a(2)))<1.e-14
                    c=polylog(a(2),exp(-1/a(1)));
                else
                    c=exp(-1/a(1));
                    inc=c;
                    j=2;
                    while inc/c>1.e-14
                        inc=j^(-a(2))*exp(-j/a(1));
                        c=c+inc;
                        j=j+1;
                    end
                end
                
                a(3)=c;
            end
            
            %now calculate the probability
            if i>=1
                y=i^(-a(2))*exp(-i/a(1))/c;
            end
            
        elseif m==1
            y=probsD(type,a,i+1,2);
            
        elseif m==2
            if length(a)>3 %if constant already calculated
                c1=a(4);   % use it
            else           %else calculate it
                if abs(a(2)-round(a(2)))<1.e-10
                    c1=polylog(a(2)-1,exp(-1/a(1)));
                else
                    c1=0*exp(-1/a(1));
                    inc=c1;
                    j=1;
                    while inc/c1>1.e-14 || j==1
                        inc=j^(-a(2)+1)*exp(-j/a(1));
                        c1=c1+inc;
                        j=j+1;
                    end
                end
                
                a(4)=c1;
            end
            
            %now calculate the probability
            if i>=1
                y=i^(-a(2)+1)*exp(-i/a(1))/c1;
            end
        end
        
        
        %My heavy-tailed thing, const for i=0,...,k then like i^(-a)
        % a(1)=k, a(2)=a
    case 'heavy'
        if m==0
            if length(a)>2 && a(3)>0 %if const already there
                c=a(3);              %  use it
            else                     %else calculate it
                c=a(1)^(-a(2)+1);
                inc=c;
                j=0;
                while inc/c>1.e-14
                    inc=(a(1)+j)^(-a(2));
                    c=c+inc;
                    j=j+1;
                end
                
                a(3)=c;
            end
            
            %now calculate the probability
            if i<=a(1)
                y=a(1)^(-a(2))/c;
            else
                y=i^(-a(2))/c;
            end
            
        elseif m==1
            y=probsD(type,a,i+1,2);
        elseif m==2
            
            if length(a)>3 %if const there
                c1=a(4);   %  use it
            else           %else calculate it
                c1=a(1)^(-a(2)+1)*(a(1)-1)/2;
                inc=c1;
                j=0;
                while inc/c1>1.e-14 || j==0
                    inc=(a(1)+j)^(-a(2)+1);
                    c1=c1+inc;
                    j=j+1;
                end
                
                a(4)=c1;
            end
            
            %calculate the probability
            if i<=a(1)-1
                y=i*a(1)^(-a(2))/c1;
            else
                y=i^(-a(2)+1)/c1;
            end
        end
        
        
        % Constant
    case 'const'
        if m==0 && i==a(1)
            y=1;
        elseif m==1 && i==a(1)-1
            y=1;
        elseif m==2 && i==a(1)
            y=1;
        end
        
        % Poisson
    case 'poi'
        if m==0 || m==1
            y=exp(-a(1))*a(1)^i/factorialbb(i);
        elseif m==2
            if i>0
                y=exp(-a(1))*a(1)^(i-1)/factorialbb(i-1);
            end
        end
        
        % Geometric
    case 'geom'
        if m==0
            y=(1-a(1))*a(1)^i;
        elseif m==1
            y=(i+1)*a(1)^i*(1-a(1))^2;
        elseif m==2
            y=i*a(1)^(i-1)*(1-a(1))^2;
        end
        
        %given
    case 'given'
        if m==1
            y=probsD(type,a,i+1,2);
            varargout{1}=y;
            return
        end
        if i>=length(a) || i<0
            y=0;
        elseif m==0
            y=a(i+1);
        elseif m==2
            y=i*a(i+1)/sum((0:(length(a)-1)).*a);
        end
        
        %Zeta
    case 'zeta'
        if m==0
            if i>0
                y=i^-(a(1))/zeta(a(1));
            end
        elseif m==2
            if i>0
                y=i^-(a(1)-1)/zeta(a(1)-1);
            end
            
        elseif m==1
            y=(i+1)^-(a(1)-1)/zeta(a(1)-1);
        end
        
        
        
        
        
        
        %Logarithmic
    case 'log'
        error('sirl:probsD:noCode','not yet')
        if m==0
            y=log(1-a(1)*s)/log(1-a(1));
        elseif m==1
            y=(1-a(1))/(1-a(1)*s);
        elseif m==2
            y=(1-a(1))/(1/s-a(1));
        end
        
        %Shifted: move a distribution to the right
    case 'shifted'
        type1=a{1};
        param1=a{2};
        shift=a{3};
        if m==0
            y=probsD(type1,param1,i-shift,0);
        elseif m==1
            y=probsD(type,a,i+1,2);
        elseif m==2
            mu=EVarD(type1,param1);
            y=i*probsD(type1,param1,i-shift,0)/(mu+shift);
        end
        
    otherwise
        error('sirl:probsD:badDistn','Unrecognised distn type')
end

varargout{1}=y;

%if the call is expecting paramD to be modified then give back the modified
%version
if nargout>1
    varargout{2}=a;
end