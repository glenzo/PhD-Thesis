#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:39:57 2021

@author: agatelucas
"""
import Bio 


dH_helixinit =  { "G" : 0.1,
                  "C" : 0.1,
                  "A" : 2.3,
                  "T" : 2.3
                }
   

dS_helixinit =  { "G" : -2.8,
                  "C" : -2.8,
                  "A" :  4.1,
                  "T" :  4.1
                }
dG_helixinit =  { "G" : 0.98,
                  "C" : 0.98,
                  "A" :  1.03,
                  "T" :  1.03
                }

dS_sym = -1.4


dH = {  "AATT" : -7.9  , "TTAA" :  -7.9,
        "ATTA" : -7.2  , "TAAT" :  -7.2,
        "CAGT" : -8.5  , "TGAC" :  -8.5,
        "GTCA" : -8.4  , "ACTG" :  -8.4,
        "CTGA" : -7.8  , "AGTC" :  -7.8,
        "GACT" : -8.2  , "TCAG" :  -8.2,
        "CGGC" : -10.6 , "GCCG" :  -9.8,
        "GGCC" : -8.0  , "CCGG" :  -8.0,
           
        # Like pair mismatches
           
        "AATA" :  1.2  , "ATAA" :   1.2,
        "CAGA" : -0.9  , "AGAC" :  -0.9,
        "GACA" : -2.9  , "ACAG" :  -2.9,
        "TAAA" :  4.7  , "AAAT" :   4.7,
       
        "ACTC" :  0.0  , "CTCA" :   0.0,
        "CCGC" : -1.5  , "CGCC" :  -1.5,
        "GCCC" :  3.6  , "CCCG" :   3.6,
        "TCAC" :  6.1  , "CACT" :   6.1,
       
        "AGTG" :-3.1   , "GTGA" :  -3.1,
        "CGGG" :-4.9   , "GGGC" :  -4.9,
        "GGCG" :-6.0   , "GCGG" :  -6.0,
        "TGAG" :1.6    , "GAGT" :   1.6,
       
        "ATTT" : -2.7  , "TTTA" : -2.7,
        "CTGT" : -5.0  , "TGTC" : -5.0,
        "GTCT" : -2.2  , "TCTG" : -2.2,
        "TTAT" : 0.2   , "TATT" :  0.2,
       
        # G.T mismatches
       
        "AGTT" :   1.0  , "TTGA" :   1.0,
        "ATTG" :  -2.5  , "GTTA" :  -2.5,
        "CGGT" :  -4.1  , "TGGC" :  -4.1,
        "CTGG" :  -2.8  , "GGTC" :  -2.8,
        "GGCT" :   3.3  , "TCGG" :   3.3,
        "GGTT" :   5.8  , "TTGG" :   5.8,
        "GTCG" :  -4.4  , "GCTG" :  -4.4,
        "GTTG" :   4.1  , "GTTG" :   4.1,
        "TGAT" :  -0.1  , "TAGT" :  -0.1,
        "TGGT" :  -1.4  , "TGGT" :  -1.4,
        "TTAG" :  -1.3  , "GATT" :  -1.3,
       
        # G.A mismatches
       
        "AATG" :  -0.6  , "GTAA" :  -0.6,
        "AGTA" :  -0.7  , "ATGA" :  -0.7,
        "CAGG" :  -0.7  , "GGAC" :  -0.7,
        "CGGA" :  -4.0  , "AGGC" :  -4.0,
        "GACG" :  -0.6  , "GCAG" :  -0.6,
        "GGCA" :   0.5  , "ACGG" :   0.5,
        "TAAG" :   0.7  , "GAAT" :   0.7,
        "TGAA" :   3.0  , "AAGT" :   3.0,
       
        # C.T mismatches
       
        "ACTT"  :  0.7  , "TTCA" :   0.7,
        "ATTC"  : -1.2  , "CTTA" :  -1.2,
        "CCGT"  : -0.8  , "TGCC" :  -0.8,
        "CTGC"  : -1.5  , "CGTC" :  -1.5,
        "GCCT"  :  2.3  , "TCCG" :   2.3,
        "GTCC"  :  5.2  , "CCTG" :   5.2,
        "TCAT"  :  1.2  , "TACT" :   1.2,
        "TTAC"  :  1.0  , "CATT" :   1.0,
       
        # A.C mismatches
       
        "AATC" :   2.3  , "CTAA" :   2.3,
        "ACTA" :   5.3  , "ATCA" :   5.3,
        "CAGC" :   1.9  , "CGAC" :   1.9,
        "CCGA" :   0.6  , "AGCC" :   0.6,
        "GACC" :   5.2  , "CCAG" :   5.2,
        "GCCA" :  -0.7  , "ACCG" :  -0.7,
        "TAAC" :   3.4  , "CAAT" :   3.4,
        "TCAA" :   7.6  , "AACT" :   7.6

}




dS = {  "AATT" : -22.2,  "TTAA" : -22.2,
        "ATTA" : -20.4,  "TAAT" : -21.3,
        "CAGT" : -22.7,  "TGAC" : -22.7,
        "GTCA" : -22.4,  "ACTG" : -22.4,
        "CTGA" : -21.0,  "AGTC" : -21.0,
        "GACT" : -22.2,  "TCAG" : -22.2,
        "CGGC" : -27.2,  "GCCG" : -24.4,
        "GGCC" : -19.9,  "CCGG" : -19.9,
       
       
        # Like pair mismatches
           
        "AATA" :   1.7,  "ATAA" :   1.7,
        "CAGA" :  -4.2,  "AGAC" :  -4.2,
        "GACA" :  -9.8,  "ACAG" :  -9.8,
        "TAAA" :  12.9,  "AAAT" :  12.9,
       
        "ACTC" :  -4.4,  "CTCA" :  -4.4,
        "CCGC" :  -7.2,  "CGCC" :  -7.2,
        "GCCC" :   8.9,  "CCCG" :   8.9,
        "TCAC" :  16.4,  "CACT" :  16.4,
       
        "AGTG" :  -9.5,  "GTGA" :  -9.5,
        "CGGG" : -15.3,  "GGGC" : -15.3,
        "GGCG" : -15.8,  "GCGG" : -15.8,
        "TGAG" :   3.6,  "GAGT" :   3.6,
       
        "ATTT" : -10.8,  "TTTA" : -10.8,
        "CTGT" : -15.8,  "TGTC" : -15.8,
        "GTCT" :  -8.4,  "TCTG" :  -8.4,
        "TTAT" :  -1.5,  "TATT" :  -1.5,
       
        # G.T mismatches
       
        "AGTT" :   0.9,  "TTGA" :   0.9,
        "ATTG" :  -8.3,  "GTTA" :  -8.3,
        "CGGT" : -11.7,  "TGGC" : -11.7,
        "CTGG" :  -8.0,  "GGTC" :  -8.0,
        "GGCT" :  10.4,  "TCGG" :  10.4,
        "GGTT" :  16.3,  "TTGG" :  16.3,
        "GTCG" : -12.3,  "GCTG" : -12.3,
        "GTTG" :   9.5,  "GTTG" :   9.5,
        "TGAT" :  -1.7,  "TAGT" :  -1.7,
        "TGGT" :  -6.2,  "TGGT" :  -6.2,
        "TTAG" :  -5.3,  "GATT" :  -5.3,
       
        # G.A mismatches
       
        "AATG" :  -2.3,  "GTAA" :  -2.3,
        "AGTA" :  -2.3,  "ATGA" :  -2.3,
        "CAGG" :  -2.3,  "GGAC" :  -2.3,
        "CGGA" : -13.2,  "AGGC" : -13.2,
        "GACG" :  -1.0,  "GCAG" :  -1.0,
        "GGCA" :   3.2,  "ACGG" :   3.2,
        "TAAG" :   0.7,  "GAAT" :   0.7,
        "TGAA" :   7.4,  "AAGT" :   7.4,
       
        # C.T mismatches
       
        "ACTT" :   0.2,  "TTCA" :   0.2,
        "ATTC" :  -6.2,  "CTTA" :  -6.2,
        "CCGT" :  -4.5,  "TGCC" :  -4.5,
        "CTGC" :  -6.1,  "CGTC" :  -6.1,
        "GCCT" :   5.4,  "TCCG" :   5.4,
        "GTCC" :  13.5,  "CCTG" :  13.5,
        "TCAT" :   0.7,  "TACT" :   0.7,
        "TTAC" :   0.7,  "CATT" :   0.7,
       
        # A.C mismatches
       
        "AATC" :   4.6,  "CTAA" :   4.6,
        "ACTA" :  14.6,  "ATCA" :  14.6,
        "CAGC" :   3.7,  "CGAC" :   3.7,
        "CCGA" :  -0.6,  "AGCC" :  -0.6,
        "GACC" :  14.2,  "CCAG" :  14.2,
        "GCCA" :  -3.8,  "ACCG" :  -3.8,
        "TAAC" :   8.0,  "CAAT" :   8.0,
        "TCAA" :  20.2,  "AACT" :  20.2
   
}
   
dH_dang = {"AA-T": 0.2, "T-AA": 0.2, "AC-G": -6.3, "G-CA": -6.3,
            "AG-C": -3.7, "C-GA": -3.7, "AT-A": -2.9, "A-TA": -2.9,
            "CA-T": 0.6, "T-AC": 0.6, "CC-G": -4.4, "G-CC": -4.4,
            "CG-C": -4.0, "C-GC": -4.0, "CT-A": -4.1, "A-TC": -4.1,
            "GA-T":-1.1, "T-AG": -1.1, "GC-G": -5.1, "G-CG": -5.1,
            "GG-C": -3.9, "C-GG": -3.9, "GT-A": -4.2, "A-TG": -4.2,
            "TA-T": -6.9, "T-AT": -6.9, "TC-G": -4.0, "G-CT": -4.0,
            "TG-C": -4.9, "C-GT": -4.9, "TT-A": -0.2, "A-TT": -0.2,
            "AAT-": -0.5, "-TAA": -0.5, "ACT-": 4.7, "-TCA": 4.7,
            "AGT-": -4.1, "-TGA": -4.1, "ATT-": -3.8, "-TTA": -3.8,
            "CAG-": -5.9, "-GAC": -5.9, "CCG-": -2.6, "-GCC": -2.6,
            "CGG-": -3.2, "-GGC": -3.2, "CTG-": -5.2, "-GTC": -5.2,
            "GAC-": -2.1, "-CAG": -2.1, "GCC-": -0.2, "-CCG": -0.2,
            "GGC-": -3.9, "-CGG": -3.9, "GTC-": -4.4, "-CTG": -4.4,
            "TAA-": -0.7, "-AAT": -0.7, "TCA-": 4.4, "-ACT": 4.4,
            "TGA-": -1.6, "-AGT": -1.6, "TTA-": 2.9, "-ATT": 2.9 }

dG_dang = { "AA-T": -0.51, "T-AA": -0.51, "AC-G": -0.96, "G-CA": -0.96,
            "AG-C": -0.58, "C-GA": -0.58, "AT-A": -0.5, "A-TA": -0.5,
            "CA-T": -0.42, "T-AC": -0.42, "CC-G": -0.52, "G-CC": -0.52,
            "CG-C": -0.34, "C-GC": -0.34, "CT-A": -0.02, "A-TC": -0.02,
            "GA-T":-0.62, "T-AG": -0.62, "GC-G": -0.72, "G-CG": -0.72,
            "GG-C": -0.56, "C-GG": -0.56, "GT-A": 0.48, "A-TG": 0.48,
            "TA-T": -0.71, "T-AT": -0.71, "TC-G": -0.58, "G-CT": -0.58,
            "TG-C": -0.61, "C-GT": -0.61, "TT-A": -0.1, "A-TT": -0.1,
            "AAT-": -0.12, "-TAA": -0.12, "ACT-": 0.28, "-TCA": 0.28,
            "AGT-": -0.01, "-TGA": -0.01, "ATT-": 0.13, "-TTA": 0.13,
            "CAG-": -0.82, "-GAC": -0.82, "CCG-": -0.31, "-GCC": -0.31,
            "CGG-": -0.01, "-GGC": -0.01, "CTG-": -0.52, "-GTC": -0.52,
            "GAC-": -0.92, "-CAG": -0.92, "GCC-": -0.23, "-CCG": -0.23,
            "GGC-": -0.44, "-CGG": -0.44, "GTC-": -0.35, "-CTG": -0.35,
            "TAA-": -0.48, "-AAT": -0.48, "TCA-": -0.19, "-ACT": -0.19,
            "TGA-": -0.5, "-AGT": -0.5, "TTA-": -0.29, "-ATT": -0.29 }

dS_dang = { "AA-T": 2.3, "T-AA": 2.3, "AC-G": -17.1, "G-CA": -17.1,
            "AG-C": -10, "C-GA": -10, "AT-A": -7.6, "A-TA": -7.6,
            "CA-T": 3.3, "T-AC": 3.3, "CC-G": -12.6, "G-CC": -12.6,
            "CG-C": -11.9, "C-GC": -11.9, "CT-A": -13, "A-TC": -13,
            "GA-T":-1.6, "T-AG": -1.6, "GC-G": -14, "G-CG": -14,
            "GG-C": -10.9, "C-GG": -10.9, "GT-A": -15, "A-TG": -15,
            "TA-T": -20, "T-AT": -20, "TC-G": -10.9, "G-CT": -10.9,
            "TG-C": -13.8, "C-GT": -13.8, "TT-A": -0.5, "A-TT": -0.5,
            "AAT-": -1.1, "-TAA": -1.1, "ACT-": 14.2, "-TCA": 14.2,
            "AGT-": -13.1, "-TGA": -13.1, "ATT-": -12.6, "-TTA": -12.6,
            "CAG-": -16.5, "-GAC": -16.5, "CCG-": -7.4, "-GCC": -7.4,
            "CGG-": -10.4, "-GGC": -10.4, "CTG-": -15, "-GTC": -15,
            "GAC-": -3.9, "-CAG": -3.9, "GCC-": -0.1, "-CCG": -0.1,
            "GGC-": -11.2, "-CGG": -11.2, "GTC-": -13.1, "-CTG": -13.1,
            "TAA-": -0.8, "-AAT": -0.8, "TCA-": 14.9, "-ACT": 14.9,
            "TGA-": -3.6, "-AGT": -3.6, "TTA-": 10.4, "-ATT": 10.4 }

internal_loop = [0,0,0,3.2,3.6,4.0,4.4,4.6,4.8,4.9,4.9]
           
R = 1.9872 # universal gas constant in Cal/degrees C*Mol
N_a = 6.022e23 #avogadro number

import numpy as np

Salt = 50e-3 # Concentration en sel (Cation monovalent, M)
SaltDiv = 1.5e-3 # Concentration en Cation divalent (M)
dNTP = 0.2e-3 # Concentration en dNTP (M)
OligoConc = 400e-9 # Concentration primer (M)

Salt = Salt+3.795*np.sqrt(SaltDiv)


if __name__ == '__main__':


        # This is just an example to play with sequences and test DG calculation


        #print(dH_helixinit)
        #print(dS_helixinit)
        #print(dH)
        #print(dS)

        import itertools

        permutations_couples = [''.join(comb) for comb in itertools.product("AGTC", repeat=4)]
        print(permutations_couples)

        for couple in permutations_couples:
                try:
                        print(dS[(couple)])
                except KeyError:
                        dS[(couple)] =0
                        dH[(couple)] =0

                       

       
from Bio import pairwise2  
       
def comp_simp(a):
   
    b = ""
   
    for n in range(len(a)):
        if a[n] == "A":
            b = b + "T"
        if a[n] == "T":
            b = b + "A"
        if a[n] == "G":
            b = b + "C"
        if a[n] == "C":
            b = b + "G"
    return b

def rev_compl(a):
    return comp_simp(a)[::-1]

def melting_temp(primer, pc1, pc2, msc, dsc, dntp, T):
   
   
    primer=primer.upper()
   
    if dsc > dntp:
        sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    else:
        sc = msc
   
    rev_primer= comp_simp(primer)
   
    dH0 = dH_helixinit[primer[0]] + dH_helixinit[primer[-1]]

    dS0 = dS_helixinit[primer[0]] + dS_helixinit[primer[-1]]
   
    k=0
   
    for p in range(len(primer)-1):      
        d_N = primer[p:p+2] + rev_primer[p:p+2]
       
        dH0 += dH[d_N]
        dS0 += dS[d_N]
        if primer[p]=="G" or primer[p] == "C" :
            k = k+1
     
     #à adapter (ne fonctionne qu'au début de la PCR)
   
    DS0=dS0
   
    """if b:
        f = k/(len(primer))
        DS0 += dS_sym/1000
        T_temp= ((dH0*1000)/(DS0 + R*np.log(c)))
    else:
        T_temp = ((dH0*1000)/(DS0 + R*np.log(c/4)))  
       
    a = 1/T_temp + (4.29*f - 3.95)*1e-5*(np.log(sc)) + 9.4e-6*(np.log(sc)**2)
    T_m2 = 1/a - 273.15 """ #température modifiée ave modèle de Van Ahsen
   
    dS0 += 0.368*(len(primer)-1)*np.log(sc)
   
    T_m = ((dH0*1000)/(dS0 + R*np.log(pc1 - pc2/2))) - 273.15  
   
       
   
       
   
    return(T_m, dH0, dS0, DS0, dH0 - (T + 273.15)*(dS0/1000), dH0 - (T + 273.15)*((dS0+ R*np.log(pc1 - pc2/2))/1000))

def melting_temp_simp(primer, pc1, pc2, msc, dsc, dntp, T):
   
   
    primer=primer.upper()
   
    if dsc > dntp:
        sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    else:
        sc = msc
   
    rev_primer= comp_simp(primer)
   
    dH0 = dH_helixinit[primer[0]] + dH_helixinit[primer[-1]]

    dS0 = dS_helixinit[primer[0]] + dS_helixinit[primer[-1]]
   
    k=0
   
    for p in range(len(primer)-1):      
        d_N = primer[p:p+2] + rev_primer[p:p+2]
       
        dH0 += dH[d_N]
        dS0 += dS[d_N]
        if primer[p]=="G" or primer[p] == "C" :
            k = k+1
     
     #à adapter (ne fonctionne qu'au début de la PCR)
   
    DS0=dS0
   
    """if b:
        f = k/(len(primer))
        DS0 += dS_sym/1000
        T_temp= ((dH0*1000)/(DS0 + R*np.log(c)))
    else:
        T_temp = ((dH0*1000)/(DS0 + R*np.log(c/4)))  
       
    a = 1/T_temp + (4.29*f - 3.95)*1e-5*(np.log(sc)) + 9.4e-6*(np.log(sc)**2)
    T_m2 = 1/a - 273.15 """ #température modifiée ave modèle de Van Ahsen
   
    dS0 += 0.368*(len(primer)-1)*np.log(sc)
   
    T_m = ((dH0*1000)/(dS0 + R*np.log(pc1 - pc2/2))) - 273.15  
   
       
   
       
   
    return(T_m, dH0 - (T + 273.15)*(dS0/1000), dH0)

def est_complementaire(a,b):
    c = False
    if a=="A" and b=="T":
        c= True
    if a=="T" and b=="A":
        c=True
    if a=="G" and b=="C":
        c=True
    if a=="C" and b=="G":
        c= True
    return c

def calcul_bulge(a,b,T):#a[0] et b[0] complémentaires, un tiret en a[1]
   
    n=len(a)
    i=1
    e=0
    s=0
   
   
    while i<n and not est_complementaire(a[i], b[i]):
        i=i+1
    i=i-1
    if i==n-1:
        e=0 #cas des dangling ends non pris en compte ici
    elif i==1:
        d_N=a[0]+a[2] + b[0] + b[2]
        dH0=dH[d_N]
        dS0=dS[d_N]
        e=4.0 + dH0 - ((T+273.15)*dS0)/1000
        s=-e*(1000/310.15)
    elif i==2:
        e=2.9
        s=-e*(1000/310.15)
    elif i==3:
        e=3.1
        s=-e*(1000/310.15)
    elif i==4:
        e=3.2
        s=-e*(1000/310.15)
    elif i==5:
        e=3.3
        s=-e*(1000/310.15)
    elif i==6:
        e=3.5
        s=-e*(1000/310.15)
    elif i==7:
        e=3.7
        s=-e*(1000/310.15)
    elif i==8:
        e=3.9
        s=-e*(1000/310.15)
    elif i==9:
        e=4.1
        s=-e*(1000/310.15)
    elif i==10:
        e=4.3
        s=-e*(1000/310.15)
    return e,i+1,0,s #i taille de la bulge loop + indice où se termine la loop

def calcul_internal(a,b,T): #les termes a[0] et b[0] sont complémentaires, a[1] et b[1] mismatch
   
    n=len(a)
    i=0
    k=1
    la=0
    lb=0
    e=0
    h=0
    s=0
    while k<(n-1) and ((not est_complementaire(a[k],b[k])) or (not est_complementaire(a[k+1],b[k+1]))):
        if a[k]=="-" or b[k] =="-":
            i=i+1
        else:
            i=i+2
        if a[k] != "-":
            la = la +1
               
        if b[k] != "-":
            lb =lb + 1
                   
        k=k+1
   
    if k>=n-1:
        e= 0
    elif i==2:
        d_N=a[0]+a[1] + b[0] + b[1]
       
        dH0=dH[d_N]
        h=h+dH0
        dS0=dS[d_N]
        s=s+dS0
        e=dH0 - ((T+273.15)*dS0)/1000
        d_N=a[1]+a[2] + b[1] + b[2]
        dH0=dH[d_N]
        h=h+dH0
        dS0=dS[d_N]
        s=s+dS0
        e= e+dH0 - ((T+273.15)*dS0)/1000
    else:
        d_N=a[0]+a[1] + b[0] + b[1]
        dH0=dH[d_N]
        dS0=dS[d_N]
        h=h+dH0
        s=s+dS0
        e=dH0 - ((T+273.15)*dS0)/1000
        d_N=a[k]+a[k+1] + b[k] + b[k+1]
        dH0=dH[d_N]
        h=h+dH0
        dS0=dS[d_N]
        s=s+dS0
        e= e+dH0 - ((T+273.15)*dS0)/1000
        if i>10:
            e = e+ dH0 - ((T+273.15)*dS0)/1000 + 0.3*abs(la-lb) + 4.9 + (10-i)*0.135
            s = s - (e*1000)/310.15
        else:
            e = e+ dH0 - ((T+273.15)*dS0)/1000 + 0.3*abs(la-lb) + internal_loop[i]
            s = s - (e*1000)/310.15
    return e,k,h,s

def complementaire(a): #complementaire mais pas dans le même sens
    n=len(a)
    b=[]
    for i in range(n):
        if a[i]=="A":
            b.append("T")
        if a[i]=="C":
            b.append("G")
        if a[i]=="G":
            b.append("C")
        if a[i]=="T":
            b.append("A")
        if a[i]=="-":
            b.append("-")
    return "".join(b)

   

def primer_dimer_simple(primer1, primer2, pc1, pc2, msc, dsc, dntp, T) : #ne considère pas dangling ends
    primer1=primer1.upper()
    primer2=primer2.upper()
    primer2=complementaire(primer2)
    a = pairwise2.align.globalms(primer1, primer2, 2, -1, -2, -.5)
    ne =len(a[0][0])
    primer1=a[0][0]
    primer2=a[0][1];primer2=complementaire(primer2)
    e=0
    h=0
    s=0
    i=0
   
    while i<(ne-1) and ((not est_complementaire(primer1[i], primer2[i])) or (est_complementaire(primer1[i], primer2[i]) and ((primer1[i+1] == "-") or (primer2[i+1] == "-")))):
        i=i+1
       
    p=i
   
    while p<(ne-1):      
        d_N = primer1[p:p+2] + primer2[p:p+2]
       
        if primer1[p+1]=="-" :
            f,d,j,t = calcul_bulge(primer1[p::], primer2[p::], T)
            p=p+d
            e=e+f
            h=h+j
            s=s+t
        elif primer2[p+1]=="-" :
            f,d,j,t = calcul_bulge(primer2[p::], primer1[p::], T)
            p=p+d
            e=e+f
            h=h+j
            s=s+t
        elif not est_complementaire(primer1[p+1], primer2[p+1]):
            f,d,j,t=calcul_internal(primer1[p::], primer2[p::],T)
            p=p+d
            e=e+f
            h=h+j
            s=s+t
        else:
            dH0 = dH[d_N]
            h=h+dH0
            dS0 = dS[d_N]
            s=s+dS0
            e=e+dH0 - ((T+273.15)*dS0)/1000
            p=p+1
       
   
    a,b,c=ajout_dang(primer1,primer2,T)
    e=e+a; h=h+b; s=s+c
   
    if dsc > dntp:
        sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    else:
        sc = msc
   
    e_60 = h - (60+273.15)*(s/1000) #((s+0.368*(len(primer1)-1)*np.log(sc) + R*np.log(pc1-pc2/2))/1000)
    #e_60 = e_60 - 0.1226*(len(primer1)-1)*np.log(sc) #- 0.662*np.log(pc1-pc2/2)
    #TODO rajouter concentration
    return primer1,primer2,e,h,s,e_60


def ajout_dang(primer1, primer2, T):
    e=0
    h=0
    s=0
    if primer1[0] == "-":
        i=0
        while i<len(primer1) and primer1[i]=="-":
            i=i+1
        if i > 0 and i<len(primer1) and est_complementaire(primer1[i], primer2[i]):
            i=i-1
            """d_N=primer1[i]+primer1[i+1] + primer2[i] + primer2[i+1]
            dH0=dH_dang[d_N]
            dG0=dG_dang[d_N]
            dS0=dS_dang[d_N]
            e=e+dG0
            h=h+dH0
            s=s+dS0"""
            h = h + dH_helixinit[primer1[i+1]]
            s = s + dS_helixinit[primer1[i+1]]
            e = e + dG_helixinit[primer1[i+1]]
    elif primer2[0] == "-":
        i=0
        while i<len(primer2) and primer2[i]=="-" :
            i=i+1
        if i > 0 and i < len(primer2) and est_complementaire(primer1[i], primer2[i]):
            i=i-1
            """d_N=primer1[i]+primer1[i+1] + primer2[i] + primer2[i+1]
            dH0=dH_dang[d_N]
            dG0=dG_dang[d_N]
            dS0=dS_dang[d_N]
            e=e+dG0
            h=h+dH0
            s=s+dS0"""
            h = h + dH_helixinit[primer1[i+1]]
            s = s + dS_helixinit[primer1[i+1]]
            e = e + dG_helixinit[primer1[i+1]]
    else:
        h = h + dH_helixinit[primer1[0]]
        s = s + dS_helixinit[primer1[0]]
        e = e + dG_helixinit[primer1[0]]
    if primer1[-1] == "-":
        i=1
        while i<len(primer1) and primer1[-i]=="-":
            i=i+1
        if i >= 0 and i<len(primer1) and est_complementaire(primer1[-i], primer2[-i]):
            i=i-1
            """d_N=primer1[-i-1]+primer1[-i] + primer2[-i-1] + primer2[-i]
            dH0=dH_dang[d_N]
            dG0=dG_dang[d_N]
            dS0=((dH0 - dG0)/(T+273.15))*1000
            e=e+dG0
            h=h+dH0
            s=s+dS0"""
            h = h + dH_helixinit[primer1[-i-1]]
            s = s + dS_helixinit[primer1[-i-1]]
            e = e + dG_helixinit[primer1[-i-1]]
    elif primer2[-1] == "-":
        i=1
        while i<len(primer2) and primer2[-i]=="-":
            i=i+1
        if i >= 0 and i<len(primer2) and est_complementaire(primer1[-i], primer2[-i]):
            i=i-1
            """d_N=primer1[-i-1]+primer1[-i] + primer2[-i-1] + primer2[-i]
            dH0=dH_dang[d_N]
            dG0=dG_dang[d_N]
            dS0=((dH0 - dG0)/(T+273.15))*1000
            e=e+dG0
            h=h+dH0
            s=s+dS0"""
            h = h + dH_helixinit[primer1[-i-1]]
            s = s + dS_helixinit[primer1[-i-1]]
            e = e + dG_helixinit[primer1[-i-1]]
    else:
        h = h + dH_helixinit[primer1[-1]]
        s = s + dS_helixinit[primer1[-1]]
        e = e + dG_helixinit[primer1[-1]]
    return e,h,s
           
def primer_dimer(primer1, primer2, pc1, pc2, msc, dsc, dntp, T):
    primer2=primer2[::-1]
   
    n=len(primer1)
    m=len(primer2)
    x = (min(n,m))*2 - 1
    l=x*[("", "", 0,0,0,0)]
    l[0]=primer_dimer_simple(primer1, primer2, pc1, pc2, msc, dsc, dntp,T)
    for i in range(1, min(n,m)):
        q=primer1[n-i::]
        r=primer2[0:i]
        l[i]=primer_dimer_simple(q,r,pc1, pc2, msc, dsc, dntp,T)
       
    for i in range (1, min(n,m)):
        q=primer2[m-i::]
        r=primer1[0:i]
        l[x-i-1]=primer_dimer_simple(q,r,pc1, pc2, msc, dsc, dntp,T)
    y=l[0][5]; t=l[0][3]; w=l[0][4]; u=l[0][2]
    z=(l[0][0], l[0][1])
    for i in range(x):
        if l[i][5] < y:
            y=l[i][5]; t=l[i][3]; w=l[i][4]; u = l[i][2]
            z=(l[i][0], l[i][1])
    return u,t,w,z, y

def primer(primer1,primer2,T):
    return primer_dimer(primer1,primer2,T), primer_dimer(primer2,primer1,T), primer_dimer(primer1,primer1,T), primer_dimer(primer2,primer2,T)

import numpy as np

def primer_temp(primer1, primer2, pc1, pc2, msc, dsc, dntp):
    a=primer_dimer(primer1,primer2,pc1, pc2, msc, dsc, dntp, 37)
    b=primer_dimer(primer2,primer1,pc1, pc2, msc, dsc, dntp, 37)
    c=primer_dimer(primer1,primer1,pc1, pc2, msc, dsc, dntp,37)    
    d=primer_dimer(primer2,primer2,pc1, pc2, msc, dsc, dntp,37)
   
    sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
   
   #à adapter avec concentration qui change;
    S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc)
    S_d = d[2] + 0.368*(len(d[3][0])-1)*np.log(sc)
    S_b = b[2] + 0.368*(len(b[3][0])-1)*np.log(sc)
    S_c = c[2] + 0.368*(len(c[3][0])-1)*np.log(sc)
   
    G_a = (a[1]*1000 - (60+273.15)*S_a)/1000
    G_b = (b[1]*1000 - (60+273.15)*S_b)/1000
    G_c = (c[1]*1000 - (60+273.15)*S_c)/1000
    G_d = (d[1]*1000 - (60+273.15)*S_d)/1000
   
    T_ma = ((a[1]*1000)/(S_a + R*np.log(pc1 - pc2/2))) - 273.15
    T_mb = ((b[1]*1000)/(S_b + R*np.log(pc1 - pc2/2))) - 273.15
    T_mc = ((c[1]*1000)/(S_c + R*np.log(pc1-pc2/2))) - 273.15
    T_md = ((d[1]*1000)/(S_d + R*np.log(pc1-pc2/2))) - 273.15
   
   
   
    return  (a,G_a, T_ma), (b,G_b,T_mb), (c,G_c,T_mc), (d,G_d,T_md)

k=1.985e-3

def self_dimer_temp60(primer, pc1, pc2, msc, dsc, dntp):
   
   
    a = primer_dimer(primer,primer,pc1, pc2, msc, dsc, dntp, 37)
    sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc) + R*np.log(pc1 - pc2/2)
   
   
   
    G_a = (a[1]*1000 - (60+273.15)*S_a)/1000
   
   
    return G_a

   
def primer_test_temp60(primer1, primer2, pc1, pc2, msc, dsc, dntp):
    a=primer_dimer(primer1,primer2,pc1, pc2, msc, dsc, dntp, 37)
    b=primer_dimer(primer2,primer1,pc1, pc2, msc, dsc, dntp, 37)
   
   
    sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
   
    if sc == 0:
        S_a = a[2]  +R*np.log(pc1 - pc2/2)
        S_b = b[2] +R*np.log(pc1 - pc2/2)
    else:
        S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc) +R*np.log(pc1 - pc2/2)
        S_b = b[2] + 0.368*(len(b[3][0])-1)*np.log(sc)+R*np.log(pc1 - pc2/2)
   
   
    G_a = (a[1]*1000 - (60+273.15)*S_a)/1000
    G_b = (b[1]*1000 - (60+273.15)*S_b)/1000
   
    m = min(G_a,G_b)
   
   
    if (a[2]==0 and a[1]==0 and a[0]==0) and (b[2]==0 and b[1]==0 and b[0]==0): #cas où les séquences ne se lient pas du tout
        m = 5*k*(60+273.15)*np.log((1 - orth_parameters.thresh_proba)/orth_parameters.thresh_proba)
       
    return m



def deltag_self_temp(primer, pc, msc, dsc, dntp, T):
    a=primer_dimer(primer,primer,37)
   
    if dsc > dntp:
        sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    else:
        sc = msc
   
    S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc) + R*np.log(pc/4)
   
   
    G_a = (a[1]*1000 - (T+273.15)*S_a)/1000
   
   
    K_a = np.exp(-G_a/(k*(T+273.15)))
   
   
    return (G_a,1-1/(1+K_a))

def self_dimer_thresh60(primer, pc, msc, dsc, dntp,x):
    a = deltag_self_temp(primer, pc, msc, dsc, dntp, 60)
    m = a[0][0]
    if m > (k*(60+273.15)*np.log((1-x)/x)):
        return True
    else:
        return False

def deltag_temp_couple(primer1, primer2, pc, msc, dsc, dntp, T):
    a=primer_dimer(primer1,primer2,37)
    b=primer_dimer(primer2,primer1,37)
   
   
    if dsc > dntp:
        sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    else:
        sc = msc
   
    S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc) + R*np.log(pc/4)
    S_b = b[2] + 0.368*(len(b[3][0])-1)*np.log(sc)+ R*np.log(pc/4)
   
   
    G_a = (a[1]*1000 - (T+273.15)*S_a)/1000
    G_b = (b[1]*1000 - (T+273.15)*S_b)/1000
   
   
    K_a = np.exp(-G_a/(k*(T+273.15)))
    K_b = np.exp(-G_b/(k*(T+273.15)))
   
   
    return (G_a,1-1/(1+K_a)), (G_b,1-1/(1+K_b))

def deltag_couple_thresh60(primer1, primer2, pc, msc, dsc, dntp, x):
    a = deltag_temp_couple(primer1, primer2, pc, msc, dsc, dntp, 60)
    m = min (a[0][0], a[1][0])
    if m > (k*(60+273.15)*np.log((1-x)/x)):
        return True
    else:
        return False



def primer_dimer_new(primer1a, primer2b, pc1, pc2, msc, dsc, dntp, T, t,u,v) :
    primer1a=primer1a.upper()
    primer2b=primer2b.upper()
    primer2b=complementaire(primer2b)
    a = pairwise2.align.globalms(primer1a, primer2b, 1, t, u, v)
    w = len(a)
    m = min(w, 2)
   
    primer1=a[0][0]
    primer2=a[0][1];primer2=complementaire(primer2)
    x = primer_dimer_simple2(primer1, primer2, pc1, pc2, msc, dsc, dntp, T)
     
    for i in range(1,m):
        p1=a[i][0]
        p2=a[i][1];p2=complementaire(p2)
        y = primer_dimer_simple2(p1, p2, pc1, pc2, msc, dsc, dntp, T)
        if y[2] < x[2]:
            x = y
           
    if len(primer1a) < len(primer2b):
        l = len(primer1a)
        n = len(primer2b)-l + 1
        for k in range(n):
            a = pairwise2.align.globalms(primer1a, primer2b[k:k+l], 1, t, u, v)
            y = primer_dimer_simple2(a[0][0], complementaire(a[0][1]), pc1, pc2, msc, dsc, dntp, T)
            if y[2] < x[2]:
                x = y
    elif len(primer1a) > len(primer2b):
        l = len(primer2b)
        n = len(primer1a)-l + 1
        for k in range(n):
            a = pairwise2.align.globalms(primer1a[k:k+l], primer2b, 1, t, u, v)
            y = primer_dimer_simple2(a[0][0], complementaire(a[0][1]), pc1, pc2, msc, dsc, dntp, T)
            if y[2] < x[2]:
                x = y
    return x
   
def primer_dimer_simple2(primer1, primer2, pc1, pc2, msc, dsc, dntp, T) : #ne considère pas dangling ends
    ne = len(primer1)
    e=0
    h=0
    s=0
    i=0
   
    while i<(ne-1) and ((not est_complementaire(primer1[i], primer2[i])) or (est_complementaire(primer1[i], primer2[i]) and ((primer1[i+1] == "-") or (primer2[i+1] == "-")))):
        i=i+1
       
    p=i
   
    while p<(ne-1):      
        d_N = primer1[p:p+2] + primer2[p:p+2]
       
        if primer1[p+1]=="-" :
            f,d,j,t = calcul_bulge(primer1[p::], primer2[p::], T)
            p=p+d
            e=e+f
            h=h+j
            s=s+t
        elif primer2[p+1]=="-" :
            f,d,j,t = calcul_bulge(primer2[p::], primer1[p::], T)
            p=p+d
            e=e+f
            h=h+j
            s=s+t
        elif not est_complementaire(primer1[p+1], primer2[p+1]):
            f,d,j,t=calcul_internal(primer1[p::], primer2[p::],T)
            p=p+d
            e=e+f
            h=h+j
            s=s+t
        else:
            dH0 = dH[d_N]
            h=h+dH0
            dS0 = dS[d_N]
            s=s+dS0
            e=e+dH0 - ((T+273.15)*dS0)/1000
            p=p+1
       
     
    a,b,c=ajout_dang(primer1,primer2,T)
    e=e+a; h=h+b; s=s+c
   
    if dsc > dntp:
        sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
    else:
        sc = msc
     #+ 1.985e-3*np.log(pc1-pc2/2)
   
    #e_37 = h - (37+273.15)*((s + 1.985e-3*np.log(pc1-pc2/2)))/1000
    if sc == 0:
        e_60 = e - (37+273.15)*1.985e-6*np.log(pc1 - pc2/2)
    else:
        e_60 = e - 0.122*(len(primer1)-1)*np.log(sc) - (37+273.15)*1.985e-6*np.log(pc1 - pc2/2)
    #TODO rajouter concentration
    return primer1,primer2,e,h,s,e_60

def primer_dimer2(primer1, primer2, pc1, pc2, msc, dsc, dntp, T,t,u,v):
    primer2=primer2[::-1]
   
    n=len(primer1)
    m=len(primer2)
    x = (min(n,m))*2 - 1
    l=x*[("", "", 0,0,0,0)]
    l[0]=primer_dimer_new(primer1, primer2, pc1, pc2, msc, dsc, dntp,T, t,u,v)
    for i in range(1, min(n,m)):
        q=primer1[n-i::]
        r=primer2[0:i]
        l[i]=primer_dimer_new(q,r,pc1, pc2, msc, dsc, dntp,T,t,u,v)
       
    for i in range (1, min(n,m)):
        q=primer2[m-i::]
        r=primer1[0:i]
        l[x-i-1]=primer_dimer_new(q,r,pc1, pc2, msc, dsc, dntp,T,t,u,v)
    y=l[0][5]; t=l[0][3]; w=l[0][4]; u=l[0][2]
    z_1, z_2 =l[0][0], l[0][1]
    for i in range(x):
        if l[i][5] < y:
            y=l[i][5]; t=l[i][3]; w=l[i][4]; u = l[i][2]
            z_1,z_2=l[i][0], l[i][1]
    return u,t,w,z_1, z_2, y


def primer_temp2(primer1, primer2, pc1, pc2, msc, dsc, dntp,t,u,v):
    a=primer_dimer2(primer1,primer2,pc1, pc2, msc, dsc, dntp, 37, t,u,v)
    b=primer_dimer2(primer2,primer1,pc1, pc2, msc, dsc, dntp, 37,t,u,v)
    #c=primer_dimer2(primer1,primer1,pc1, pc2, msc, dsc, dntp,37,t,u,v)    
    #d=primer_dimer2(primer2,primer2,pc1, pc2, msc, dsc, dntp,37,t,u,v)
   
    sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
   
   
   #à adapter avec concentration qui change;
   
   
    S_a = a[2] + 0.368*(len(a[3])-1)*np.log(sc) #+ R*np.log(pc1 - pc2/2)
    S_b = b[2] + 0.368*(len(b[3])-1)*np.log(sc) #+ R*np.log(pc1 - pc2/2)
    #S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc)
    #S_d = d[2] + 0.368*(len(d[3][0])-1)*np.log(sc)
    #S_b = b[2] + 0.368*(len(b[3][0])-*np.log(sc)
    #S_c = c[2] + 0.368*(len(c[3][0])-1)*np.log(sc)
   
    G_a = (a[1]*1000 - (60+273.15)*(S_a))/1000
    G_b = (b[1]*1000 - (60+273.15)*(S_b))/1000
    #G_c = (c[1]*1000 - (60+273.15)*(S_c + R*np.log(pc1 - pc2/2)))/1000
    #G_d = (d[1]*1000 - (60+273.15)*(S_d + R*np.log(pc1 - pc2/2)))/1000
   
    m = min(G_a, G_b) # G_c, G_d)
   
    T_ma = ((a[1]*1000)/(S_a+R*np.log(pc1 - pc2/2))) - 273.15
    T_mb = ((b[1]*1000)/(S_b+R*np.log(pc1 - pc2/2))) - 273.15
    #T_mc = ((c[1]*1000)/(S_c + R*np.log(pc1-pc2/2))) - 273.15
    #T_md = ((d[1]*1000)/(S_d + R*np.log(pc1-pc2/2))) - 273.15
   
    x = max(T_ma, T_mb)
   
    return m, x, (a,S_a,G_a,T_ma), (b,S_b,G_b,T_mb) #,(c,G_c,T_mc), (d,G_d,T_md)

def liste_mutation(a):
    n = len(a)
    l = []
    for i in range(n):
        p1 = a
        pr1 = list(p1)
        if p1[i] == "A":
            pr1[i] = "G"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "C"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "T"
           
            pri1 = "". join(pr1)
            l.append(pri1)
        elif p1[i] == "C":
            pr1[i] = "G"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "A"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "T"
            pri1 = "". join(pr1)
            l.append(pri1)
        elif p1[i] == "T":
            pr1[i] = "G"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "C"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "A"
            pri1 = "". join(pr1)
            l.append(pri1)
        elif p1[i] == "G":
            pr1[i] = "A"
            pri1 = "". join(pr1)
            l.append(pri1)
            pr1[i] = "C"
            pri1 = "". join(pr1)
            l.append(pri1)
           
            pr1[i] = "T"
            pri1 = "". join(pr1)
            l.append(pri1)
    return l
   
   
def mutation1(primer1, primer2, pc1, pc2, msc, dsc, dntp):
    l = liste_mutation(primer1)
    n = len(l)
    s=0
    c = primer_temp2(primer1, primer2, pc1, pc2, msc, dsc, dntp, 10)[0][2]
    for y in l:
        x = primer_temp2(y, primer2, pc1, pc2, msc, dsc, dntp, 10)[0][2]
        s = s + (c-x)
    return (s/(n))

def mutation2(primer1, primer2, pc1, pc2, msc, dsc, dntp):
    l = liste_mutation(primer1)
    n = len(l)
    s=0
    c = primer_temp2(primer1, primer2, pc1, pc2, msc, dsc, dntp, 10)[0][2]
    for i in range(n):
        li = liste_mutation(l[i])
        for y in li:
            if not y == primer1:
                x = primer_temp2(y, primer2, pc1, pc2, msc, dsc, dntp, 10)[0][2]
                s = s + (c-x)
    return ((2*s)/(n*(n-1)))

def primer_temp_test(primer1, primer2, pc1, pc2, msc, dsc, dntp,t,u,v):
    a=primer_dimer2(primer1,primer2,pc1, pc2, msc, dsc, dntp, 37,t,u,v)
    b=primer_dimer2(primer2,primer1,pc1, pc2, msc, dsc, dntp, 37,t,u,v)
   
   
    sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
   
   #à adapter avec concentration qui change;
    S_a = a[2] + 0.368*(len(a[3][0])-1)*np.log(sc)
    S_b = b[2] + 0.368*(len(b[3][0])-1)*np.log(sc)
   
   
    G_a = (a[1]*1000 - (37+273.15)*(S_a+ R*np.log(pc1 - pc2/2)))/1000
    G_b = (b[1]*1000 - (37+273.15)*(S_b+ R*np.log(pc1 - pc2/2)))/1000
   

    if a[2]==0 and a[1]==0 and a[0]==0:
        G_a = 20; G_b = 20
    if G_a < G_b:
        return G_a
    else:
        return G_b
   
def primer_test60(primer1, primer2, pc1, pc2, msc, dsc, dntp,t,u,v):
    a=primer_dimer2(primer1,primer2,pc1, pc2, msc, dsc, dntp, 37,t,u,v)
    b=primer_dimer2(primer2,primer1,pc1, pc2, msc, dsc, dntp, 37,t,u,v)
    print(a)
   
   
    sc = msc + (120*np.sqrt(1000*dsc - 1000*dntp))/1000
   
   #à adapter avec concentration qui change;
    S_a = a[2] + 0.368*(len(a[3])-1)*np.log(sc)
    S_b = b[2] + 0.368*(len(b[3])-1)*np.log(sc)
   
   
    G_a = (a[1]*1000 - (60+273.15)*(S_a+ R*np.log(pc1 - pc2/2)))/1000
    G_b = (b[1]*1000 - (60+273.15)*(S_b+ R*np.log(pc1 - pc2/2)))/1000
   
   

    if a[2]==0 and a[1]==0 and a[0]==0 and b[2]==0 and b[1]==0 and b[0]==0:
        G_a = 20; G_b = 20
   
    return min(G_a, G_b)
   

def genere_sequences(n):
    s = ""
    for i in range(n):
        a = np.random.randint(1,5)
        if a == 1:
            s = s + "A"
        elif a ==2 :
            s = s + "C"
        elif a ==3:
            s = s + "G"
        elif a == 4:
            s = s+ "T"
    return s

import time

def calcul_parametre(p,n):
    li = []
    dico_tot = {}
    dico_moy = {}
    dico_rar = {}
    dico_un = {}
    t1 = time.perf_counter()
    for i in range(p):
        x = 20
        y = (0,0,0)
        seq1 = genere_sequences(n)
        seq2 = genere_sequences(n)
        lis = []
        for j in range (5):
            a = - j/2
           
            for k in range(5):
                b = -1 - k/2
               
                l = 0
                c= 0
                while (c > b) :
                    c =  - (l/2)
                    l = l + 1
                   
                    z = primer_temp_test(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, a,b,c)
                    if z <= x:
                        x = z
                        y = (a,b,c)
                        lis.append((x,y))
   
        test = sorted(lis)
       
        val = test[0][0]
        while test[-1][0] != val:
            test.pop()
       
        length = len(test)
        for (v, tri) in test:
            if (tri in dico_tot):
                dico_tot[tri] += 1
                dico_moy[tri] += 1/length
                if length <= 3 :
                    if (tri in dico_rar):
                        dico_rar[tri] +=1
                    else:
                        dico_rar[tri] = 1
                if length == 1:
                    if (tri in dico_un):
                        dico_un[tri] += 1
                    else :
                        dico_un[tri] = 1
            else:
                dico_tot[tri] = 1
                dico_moy[tri] = 1/length
                if length<=3:
                    dico_rar[tri]=1
                if length == 1 :
                    dico_un[tri] = 1
        t2 = time.perf_counter()
        print(t2 - t1)
       
    return dico_tot, dico_moy, dico_rar, dico_un

def find_primer(seq1,b):
    n = len(b)
    i = 0
    li = []
    t_pf = melting_temp(seq1, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 60)[4]
    while i + 18 < n :
        j = 18
        while (j < 25) and (i + j < n) :
            seq2 = b[i:i+j]
            test_1 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -0.5,-1,-1)
            test_2 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,-1)
            test_3 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,0)
            test_4 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-1,0)
            test_5 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -0.5,-3,0)
            test_6 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-1,-1)
            test_7 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-2,0)
            test_8 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-2,-1)
            test_9 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-2,0)
            test_10 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-2,-1)
            test_11 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-2,-1)
            test_12 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-1,-1)
            test_13 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-3,-1)
            test_14 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-2,0)
            mi = min(test_1, test_2, test_3, test_4, test_5, test_6, test_7,test_8, test_9, test_10, test_11, test_12, test_13, test_14)
            t_pr = melting_temp(seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 60)[4]
            li.append((mi, (t_pr-t_pf),seq2))
           
            j = j +1
        i = i+1    
    return sorted(li)

liste_coord = [(-0.5,-1,-1),(-2,-1,-1),(-2,-1,0),(0,-1,0),(-0.5,-3,0),(0,-1,-1),(-2,-2,0),(-2,-2,-1),(-1,-2,0),(0,-2,-1),(-1,-2,-1),(-1,-1,-1),(-1,-3,-1),(0,-2,0)]

def primer_temp_check(seq1, seq2, pc1, pc2, msc, dsc, dntp) :
       test_1 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -0.5,-1,-1)[0]
       test_2 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,-1)[0]
       test_3 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,0)[0]
       test_4 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-1,0)[0]
       test_5 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -0.5,-3,0)[0]
       test_6 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-1,-1)[0]
       test_7 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-2,0)[0]
       test_8 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-2,-1)[0]
       test_9 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-2,0)[0]
       test_10 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-2,-1)[0]
       test_11 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-2,-1)[0]
       test_12 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-1,-1)[0]
       test_13 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -1,-3,-1)[0]
       test_14 = primer_temp2(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-2,0)[0]
       li = [test_1, test_2, test_3, test_4, test_5, test_6, test_7,test_8, test_9, test_10, test_11, test_12, test_13, test_14]
       ind = li.index(min(li))
       t,u,v = liste_coord[ind]
       
       
       return primer_temp2(seq1,seq2,pc1,pc2,msc,dsc,dntp,t,u,v)
   
def primer_temp_fastcheck(seq1, seq2, pc1, pc2, msc, dsc, dntp) :
       test_1 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -0.5,-1,-1)
       test_2 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,-1)
       test_3 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,0)
       test_4 = primer_test60(seq1, seq2, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-1,0)
       
       li = [test_1, test_2, test_3, test_4]
       
       
       return min(li)
   
def primer_fast_self_check(seq1, pc1, pc2, msc, dsc, dntp) :
    test_1 = primer_test60(seq1, seq1, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -0.5,-1,-1)
    test_2 = primer_test60(seq1, seq1, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,-1)
    test_3 = primer_test60(seq1, seq1, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, -2,-1,0)
    test_4 = primer_test60(seq1, seq1, 400e-9, 400e-9, 50e-3, 3e-3, 200e-6, 0,-1,0)
       
    li = [test_1, test_2, test_3, test_4]
   
       
    return min(li)

def equivalence(msc):
   
   
    x = ((100/12)**2) * (msc**2)* (1/1000)
    return x
       
       
if __name__ == '__main__':
    print(2)
    with open('/Users/glencarter/Desktop/PhD/Scripts/thermo_codes/load_sequences/stem_blocker_seqs.txt','r') as input_file:
        with open('/Users/glencarter/Desktop//PhD/Scripts/thermo_codes/load_sequences/stem_blocker_results_file.txt','w') as output_file:
            k = 0
            for l in input_file:
                
                x = l.replace('\n', '')
                if x != '':
                    k = k+1
                    a, b, c = melting_temp_simp(x, 10e-9, 0.25e-9, 150e-3, 0, 0, 25)
                    output_file.write(x + 'T_m' + str(a) + ' DeltaG(25): ' + str(b) + '\n')
            print(k)
        """output_file.write("Dictionnaire unique" + '\n')
        for k, v in sorted(d.items(), key=lambda x: x[1]):
            output_file.write("%s: %s" % (k, v) + '\n') """

