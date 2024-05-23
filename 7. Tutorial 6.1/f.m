function fout = f(S,r_0,r_max,sig,x);

if S < 0
fout = r_0;
else
fout = r_0 + r_max*(S^x)/(S^x+sig^x );
end