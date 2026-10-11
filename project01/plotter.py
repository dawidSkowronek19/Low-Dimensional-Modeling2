import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import sys

U = []
U0 = []
S=[]
D=[]


os.makedirs('graphs', exist_ok=True)

print('n | u(0) | S | D')
for n in range(2, 40, 2):
    result = subprocess.run([sys.argv[1], str(n)], capture_output=True, text=True)
    result = result.stdout.strip().split(' | ')
    data = np.loadtxt(f'./tmp/u{n}.txt')
    S.append(float(result[2]))
    U.append(data[:,1])
    U0.append(data[np.abs(data[:,0]) < 1e-8 ,1][0])
    plt.figure(figsize=(10,8))
    plt.grid()
    plt.plot(data[:,0], data[:,1])
    plt.savefig(f'graphs/u{n}.png')
    plt.close()



    if (n<=2):
        D.append(0.0)
        continue
    else:
        D.append(np.sqrt(((U[-1]-U[-2])**2).sum()))
        print()
    print(f'{result} | {D[-1]}')


n=np.arange(2, 40, 2)
fig, axs = plt.subplots(4, 1, figsize=(14, 10))
axs[0].plot(n, U0, color='blue')
axs[0].set_title(r'u(n,x=0)')


axs[1].plot(n, S, color='red')
axs[1].set_title('S(n)')

axs[2].plot(n[1:], D[1:], color='green')
axs[2].set_title('D(n)')

x = np.linspace(-1, 1, len(U[0]))
for n_idx in range(2, 40, 8):
    axs[3].plot(x, U[(n_idx // 2) - 1], label=f'n={n_idx}')

axs[3].set_title('u(n, x)')
axs[3].legend(fontsize='small', ncol=4)
for ax in axs:
    ax.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig('./graphs/output', dpi=300)