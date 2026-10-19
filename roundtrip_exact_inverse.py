import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

a_env = 2.0
D = 0.5
a0 = a_env - 1j*D
a_D = -1j*D
alpha = np.sqrt(np.pi/a_D)
Dtilde = -1/(4*D)
N = 256

def make_grid(N, L):
    dx = L/N
    x = (np.arange(N)-N//2)*dx
    kx = np.fft.fftshift(np.fft.fftfreq(N, d=dx))*2*np.pi
    return x, dx, kx

def aspw_roundtrip(N, L):
    x, dx, kx = make_grid(N, L)
    V0 = np.exp(-a0*x**2)
    V_hat = np.fft.fftshift(np.fft.fft(np.fft.fftshift(V0)))*dx
    V_back = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(V_hat)))/dx
    return x, V0, V_back, kx, V_hat

def semi_forward_and_back(N, L):
    x, dx, kx = make_grid(N, L)
    dk = kx[1]-kx[0]
    beta = -2*Dtilde*kx
    M = np.exp(1j*np.outer(beta, kx)) * dk/(2*np.pi)     # forward "second-transform" matrix: h = M @ g
    Minv = np.linalg.inv(M)                                # its EXACT algebraic inverse

    U = np.exp(-a_env*x**2)
    U_hat = np.fft.fftshift(np.fft.fft(np.fft.fftshift(U)))*dx
    g = U_hat * np.exp(1j*Dtilde*kx**2)
    h = M @ g
    V_hat = alpha * np.exp(1j*Dtilde*kx**2) * h              # forward semi-analytical spectrum

    # ---- inverse: exact mirror using Minv ----
    A_tilde = V_hat * np.exp(-1j*Dtilde*kx**2)                # = alpha * h
    h_rec = A_tilde/alpha
    g_rec = Minv @ h_rec
    U_hat_rec = g_rec * np.exp(-1j*Dtilde*kx**2)
    U_rec = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(U_hat_rec)))/dx
    V_back = U_rec * np.exp(1j*D*x**2)

    V0 = np.exp(-a_env*x**2) * np.exp(1j*D*x**2)
    return x, V0, V_back, kx, V_hat

Ls = [1, 4, 10, 20]
fig, axes = plt.subplots(4, 2, figsize=(13, 16))

print(f"{'L':>5} {'ASPW rt err':>14} {'Semi rt err':>14} {'Semi fwd err (vs exact)':>26}")
for row, L in enumerate(Ls):
    x_a, V0_a, Vback_a, kx_a, Vhat_a = aspw_roundtrip(N, L)
    x_s, V0_s, Vback_s, kx_s, Vhat_s = semi_forward_and_back(N, L)

    Vhat_exact = np.sqrt(np.pi/a0)*np.exp(-kx_s**2/(4*a0))
    err_fwd_semi = np.linalg.norm(Vhat_s-Vhat_exact)/np.linalg.norm(Vhat_exact)

    err_a = np.linalg.norm(Vback_a-V0_a)/np.linalg.norm(V0_a)
    err_s = np.linalg.norm(Vback_s-V0_s)/np.linalg.norm(V0_s)
    print(f"{L:5.0f} {err_a:14.4e} {err_s:14.4e} {err_fwd_semi:26.4e}")

    mask = np.abs(x_a) <= min(3, L/2)
    ax_amp = axes[row,0]
    ax_amp.plot(x_a[mask], np.abs(V0_a[mask]), color='black', lw=2, label='original')
    ax_amp.plot(x_a[mask], np.abs(Vback_a[mask]), 'o--', color='crimson', ms=3, label='ASPW round-trip')
    ax_amp.plot(x_s[mask], np.abs(Vback_s[mask]), 's--', color='teal', ms=3, label='Semi-analytical round-trip')
    ax_amp.set_title(f'L={L}mm: |V(x)| round-trip')
    ax_amp.set_xlabel('x (mm)'); ax_amp.legend(fontsize=8)

    ax_ph = axes[row,1]
    ax_ph.plot(x_a[mask], np.unwrap(np.angle(V0_a[mask])), color='black', lw=2, label='original')
    ax_ph.plot(x_a[mask], np.unwrap(np.angle(Vback_a[mask])), 'o--', color='crimson', ms=3, label='ASPW round-trip')
    ax_ph.plot(x_s[mask], np.unwrap(np.angle(Vback_s[mask])), 's--', color='teal', ms=3, label='Semi-analytical round-trip')
    ax_ph.set_title(f'L={L}mm: phase(V(x)) round-trip')
    ax_ph.set_xlabel('x (mm)'); ax_ph.legend(fontsize=8)

plt.suptitle(r'Round trip (forward+inverse), $\exp(-2x^2)\exp(ix^2/2)$, N=256', y=1.0)
plt.tight_layout()
plt.savefig('/home/claude/roundtrip_exact_inverse.png', dpi=130)
