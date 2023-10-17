from dNdS import reportCI,CfracdNdS

k=30
nt_cont=0.403531
prot_cont=0.750943
PdS = CfracdNdS.calc_PdS(prot_cont,nt_cont,k)
PdN = CfracdNdS.calc_PdN(prot_cont,k)
PdS_const = PdS/0.77
PdN_const = PdN/2.23
print(PdN, PdS, PdN/PdS)
print(PdN_const, PdS_const, PdN_const/PdS_const)
