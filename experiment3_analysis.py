import json
import matplotlib.pyplot as plt
import numpy as np

#This font-set has the right kind of epsilon
plt.rcParams["mathtext.fontset"] = "cm";
plt.rcParams["font.size"] = 14;
plt.rcParams["font.weight"] = 32;

names = ['transport_optimal_BEC.json','transport_optimal_BEC_slow.json','transport_optimal_BEC_sloow.json','shake_optimal_BEC_N0.json']

algorithms= ['GRAPE_BFGS_L2','GRAPE_STEEPEST_L2','GRAPE_BFGS_H1','GRAPE_STEEPEST_H1']
algorithms_printable= ['Bfgs/L2','Steep/L2','Bfgs/H1','Steep/H1']

coeffs = json.load(open('coefficients.json'))
for name in names:
    data = json.load(open(name))

    xs = data["x"]
    ts = data["t"]

    Ns  = data["N"]

    if not isinstance(Ns,list):
        Ns = [Ns]

    for (i,N) in enumerate(Ns):
        for (j,A) in enumerate(algorithms):
            fig, ax = plt.subplots()
            psi = data['psi2_'+A+'_' + str(i)]

            ax.imshow(np.abs(psi))
            xticks = np.linspace(0,len(ts),round(ts[-1])*2+1)
            ax.set_xticks(xticks)

            if N == 0:
                ax.set_title('g=0 | '+algorithms_printable[j]);
            else:
                ax.set_title('N='+str(N)+" | "+algorithms_printable[j]);

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
        fmts = ['-.',':','-.',':']
        for (j,A) in enumerate(algorithms):
            if N == 0:
                plt.plot(ts,data['u_'+A+'_' + str(i)][0],fmts[j],label=('g=0 | '+algorithms_printable[j]))
            else:
                plt.plot(ts,data['u_'+A+'_' + str(i)][0],fmts[j],label=('N='+str(round(N))+" | "+algorithms_printable[j]))

        plt.grid()
        plt.legend()
        plt.xlabel(r'$t/ms$')
        plt.ylabel(r'$x_0/\mu s$')
        plt.show()

        fig, ax = plt.subplots()
        xticks = [1,2,3,4]
        for (j,A) in enumerate(algorithms):
            if (data['use_'+A]==1):
                plt.plot(xticks[j],data['F_'+A+'_' + str(i)],'ro')
        ax.set_xticks(xticks)
        ax.set_xlim([0,5]);
        ax.set_ylim([0,1.1]);

        ax.set_xticklabels(algorithms_printable)

        ax.grid();

        ax.set_ylabel(r'$\mathcal{F}$')
        plt.show()

        fig, ax = plt.subplots()
        xticks = [1,2,3,4]
        for (j,A) in enumerate(algorithms):
            if (data['use_'+A]==1):
                plt.plot(xticks[j],data['T_compute_'+A+'_' + str(i)],'ro')
        ax.set_xticks(xticks)
        ax.set_xlim([0,5]);

        ax.set_xticklabels(algorithms_printable)

        ax.grid();
        ax.set_ylabel(r'$T_{cacl}/s$')
        plt.show()

        fig, ax = plt.subplots()
        xticks = [1,2,3,4]
        for (j,A) in enumerate(algorithms):

            if (data['use_'+A]==1):
                plt.plot(xticks[j],data['Nsteps_'+A+'_' + str(i)],'ro')
        ax.set_xticks(xticks)
        ax.set_xlim([0,5]);

        ax.set_xticklabels(algorithms_printable)

        ax.grid();
        ax.set_ylabel('iterations')
        plt.show()
