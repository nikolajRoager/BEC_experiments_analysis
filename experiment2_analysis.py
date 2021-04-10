import json
import matplotlib.pyplot as plt
import numpy as np

#This font-set has the right kind of epsilon
plt.rcParams["mathtext.fontset"] = "cm";
plt.rcParams["font.size"] = 14;
plt.rcParams["font.weight"] = 32;

names = ['squeeze_BEC.json','shake_BEC.json', 'shake_BEC_focus.json']

coeffs = json.load(open('coefficients.json'))
for name in names:
    data = json.load(open(name))

    xs = data["x"]
    ts = data["t"]

    Ns  = data["N"]

    if not isinstance(Ns,list):
        Ns = [Ns]

    for (i,N) in enumerate(Ns):
        print(N)
        plt.plot(ts,data["Z0_"+str(i)],label=r'$Z_0$');
        plt.plot(ts,data["Z1_"+str(i)],label=r'$Z_1$');
        plt.plot(ts,data["Z2_"+str(i)],label=r'$Z_2$');
        plt.xlabel(r'$t/ms$')
        plt.ylabel(r'$Z_i$')
        plt.grid();
        plt.legend();
        plt.show();
        fig, ax = plt.subplots()


        psi = data['psi2_' + str(i)]

        ax.imshow(np.abs(psi))
        xticks = np.linspace(0,len(ts),round(ts[-1])*2+1)
        ax.set_xticks(xticks)

        ticks =[];
        for t in np.linspace(0,ts[-1],len(xticks)):
            ticks.append(str(round(t,1)));
        ax.set_xticklabels(ticks)

        yticks = np.linspace(0,len(xs),round(xs[-1])+1)
        ax.set_yticks(yticks)
        ticks =[];
        for x in np.linspace(0,xs[-1],len(yticks)):
            ticks.append(str(round(x*2-len(yticks)+1,1)));
        ax.set_yticklabels(ticks)

        ax.set_xlabel(r'$t/ms$')
        ax.set_ylabel(r'$x/\mu m$')
        plt.show()

    for (i,N) in enumerate(Ns):
        if N==0:
            plt.plot(ts,data["Z0_"+str(i)],label=r'$g=0$');
        elif (N!=10000):
            plt.plot(ts,data["Z0_"+str(i)],label=r'$N='+str(round(N))+'$');
        plt.xlabel(r'$t/ms$')
        plt.ylabel(r'$Z_0$')
    plt.grid();
    plt.legend();
    plt.show();

