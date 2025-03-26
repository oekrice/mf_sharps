import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

matplotlib.rcParams['text.usetex'] = True

fig_width = 513.11743/72

root1 = './hdata/'  #For new tests

cs = ['red', 'blue', 'green']

fig,axs = plt.subplots(2, figsize = (fig_width, fig_width*0.7))
runs = [31,32,33,34,35,36,37]

cs = plt.cm.plasma(np.linspace(0.1,0.9,len(runs)))

for ri, run in enumerate(runs):
    try:
        h_all = np.load(root1 + 'halls%03d.npy' % run)
        omegas = np.load(root1 + 'omegas%03d.npy' % run)
        tplots = np.load(root1 + 'tplots%03d.npy' % run)
        end = np.sum(np.abs(h_all) > 0) - 1
        axs[0].plot(tplots[:end], h_all[:end])
        axs[1].plot(tplots[:end], omegas[:end])
    except:
        pass

    print(h_all)
t_ref = np.load(root1 + 't_ref.npy')
h_ref = np.load(root1 + 'h_ref.npy')


axs[0].plot(t_ref, h_ref, c = 'black', linestyle = 'dashed', label = 'MHD Reference')
axs[0].set_xlim(0,t_ref[-1])
axs[0].set_title('In-plane Helicity $H$')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('$H$')

axs[1].set_xlim(0,t_ref[-1])
axs[1].set_title('Ideal Omega')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('$\Omega$')

plt.tight_layout()
plt.savefig('helicity_adapt_sharp.png')
plt.show()
