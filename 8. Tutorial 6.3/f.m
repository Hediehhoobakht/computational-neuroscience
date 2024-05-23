function fout = fin(V_ss)
Vth=-50e-3;
Vreset=-80e-3;%-0.08
sigv=1e-3;
tau=3e-3;
if V_ss==Vth;
fout=sigv/(tau*(V_ss-Vreset));
else
fout=(V_ss-Vth)/(tau*(Vth-Vreset)*(1-exp(-(V_ss-Vth)/sigv)));
end