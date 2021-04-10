import json
import matplotlib.pyplot as plt
import numpy as np

#This font-set has the right kind of epsilon
plt.rcParams["mathtext.fontset"] = "cm";
plt.rcParams["font.size"] = 14;
plt.rcParams["font.weight"] = 32;

names = ['someNs.json','manyNs.json']

coeffs = json.load(open('coefficients.json'))

A     = coeffs["A"];
sigma = coeffs["sigma"];

for name in names:

    data = json.load(open(name))


    xs = data["x"]
    V  = data["V"]


    fig,axes = plt.subplots(1,1);
    ax=axes;
    ax.grid();
    ax.set_xlabel(r'$x/\mu m$')
    ax.set_ylabel(r'$V/\epsilon$')
    ax.plot(xs,V,'r:');

    ax2=ax.twinx()
    ax2.set_ylabel(r'$|\psi|^2/\xi^2$')
    Ns  = data["N"]

    for (i,N) in enumerate(Ns):
        psi2 = data["psi2_"+str(i)]
        if N == 0:
            ax2.plot(xs,psi2,label="g=0");
        else:
            ax2.plot(xs,psi2,label="N="+str(N));

    plt.legend();
    plt.show();

    plt.grid();
    plt.xlabel(r'$N$')
    plt.ylabel(r'$\sigma^2_{|\psi|}/\chi^2$')
    plt.plot(Ns,(np.array(data["sig2"])),'r:.');

    plt.show();

    plt.grid();
    plt.xlabel(r'$N$')
    plt.ylabel(r'$T_{calc}/s$')
    plt.plot(Ns,np.sqrt(np.array(data["Tcompute"])),'r.');
    plt.show();

    #plt.grid();

    #plt.xlabel(r'$x/\mu m$')
    #plt.ylabel(r'$y/\epsilon$')

    plt.show();
