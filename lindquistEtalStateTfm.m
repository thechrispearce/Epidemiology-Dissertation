function f=lindquistEtalStateTfm(M,x,i,j)

%transform (x,i,j) into a linear index for the Lindquist etal SIS ODEs

%M is the max degree being considered
%x = 0/1 means S/I
%i,j are subscripts, non-neg integers with i+j \leq M

f = j+1 + i*(2*M+3-i)/2;

if x==1
    f=f+(M+1)*(M+2)/2;
end