import numpy as np
import matplotlib.pyplot as plt
import scipy

def angular_spectrum_method(E, z, k2, Nx, dx):
    """
    Compute the field in distance equal to z with the angular spectrum method. 
    By default (scale_factor = 1), the ouplut plane coordinates is the same than the input.
    Otherwise, it's recommended to use the two_steps_fresnel_method as it's computationally cheaper.
    To arbitrarily choose and zoom in a region of interest, use bluestein method instead.

    Reference: https://rafael-fuente.github.io/simulating-diffraction-patterns-with-the-angular-spectrum-method-and-python.html
    """
    # compute angular spectrum
    fft_c = scipy.fft.fft2(E)
    c = scipy.fft.fftshift(fft_c)

    fx = 2 * np.pi * scipy.fft.fftshift(scipy.fft.fftfreq(Nx, d = dx))
    fy = 2 * np.pi * scipy.fft.fftshift(scipy.fft.fftfreq(Nx, d = dx))
    kx, ky = np.meshgrid(fx, fy)

    argument = k2 - kx*kx - ky*ky

    #Calculate the propagating and the evanescent (complex) modes
    tmp = np.sqrt(argument)
    # kz = np.where(argument >= 0, tmp, 1j*tmp)
    kz = tmp

    # propagate the angular spectrum a distance z
    E = scipy.fft.ifft2(scipy.fft.ifftshift(c * np.exp(1j * kz* z)))
    # E = np.fft.ifft2(c)

    return E

def quadratic_taylor_terms(func, Nx, dx):
    Cxx = (func[Nx//2+1, Nx//2] - 2*func[Nx//2, Nx//2] + func[Nx//2-1, Nx//2]) / (2*dx*dx)
    Cyy = (func[Nx//2, Nx//2+1] - 2*func[Nx//2, Nx//2] + func[Nx//2, Nx//2-1]) / (2*dx*dx)
    Cxy = (func[Nx//2+1, Nx//2+1] - func[Nx//2+1, Nx//2-1] - func[Nx//2-1, Nx//2+1] + func[Nx//2-1, Nx//2-1]) / (4*dx*dx)
    return Cxx, Cxy, Cyy


def semi_analytical_asm(E, z, k2, Nx, dx, xx, yy):
    """
    Angular spectruum method from: https://rafael-fuente.github.io/simulating-diffraction-patterns-with-the-angular-spectrum-method-and-python.html

    Modified by Karthik Sreekumar
    """

    # ---------------------------
    # Forward semi analytical FFT
    # ---------------------------

    # computing quadratic and residual phase
    phase = np.atan2(E.imag, E.real)
    Cxx, Cxy, Cyy = quadratic_taylor_terms(phase, Nx, dx)
    phase_quad = Cxx*xx*xx + Cyy*yy*yy + Cxy*xx*yy
    residual_phase = phase - phase_quad

    U = np.abs(E) * np.exp(1j*residual_phase)

    # calculating analytical part of fourier transform
    fx = 2 * np.pi * scipy.fft.fftshift(scipy.fft.fftfreq(Nx, d = dx))
    fy = 2 * np.pi * scipy.fft.fftshift(scipy.fft.fftfreq(Nx, d = dx))
    kx, ky = np.meshgrid(fx, fy)

    denom = Cxy*Cxy - 4*Cxx*Cyy
    Cxx_inv = Cyy/denom
    Cyy_inv = Cxx/denom
    Cxy_inv = -Cxy/denom
    phase_quad_inv = Cxx_inv*kx*kx + Cyy_inv*ky*ky + Cxy_inv*kx*ky
    alpha = 0
    if Cxx != 0:
        alpha = 1j*np.sqrt(-1/denom + 0j)
    else:
        alpha = 1/Cxy

    # compute forward semi analytical FFT
    fft_c = scipy.fft.fft2(U)
    c = scipy.fft.fftshift(fft_c)
    U = c * np.exp(1j * phase_quad_inv)
    E = alpha * scipy.fft.ifft2(scipy.fft.ifftshift(U)) * np.exp(1j * phase_quad_inv)

    plt.imshow(np.abs(E))
    plt.colorbar()
    plt.show()

    argument = k2 - kx*kx - ky*ky

    #Calculate the propagating and the evanescent (complex) modes
    tmp = np.sqrt(argument)
    # kz = np.where(argument >= 0, tmp, 1j*tmp)
    kz = tmp

    # propagate the angular spectrum a distance z
    E = E * np.exp(1j * kz* z)
    # E = np.fft.ifft2(c)

    # ---------------------------
    # Inverse semi analytical FFT
    # ---------------------------

    # computing quadratic and residual phase
    phase = np.atan2(E.imag, E.real)
    Cxx_inv, Cxy_inv, Cyy_inv = quadratic_taylor_terms(phase, Nx, np.abs(fx[1] - fx[0]))
    phase_quad_inv = Cxx_inv*kx*kx + Cyy_inv*ky*ky + Cxy_inv*kx*ky
    residual_phase = phase - phase_quad_inv

    U = np.abs(E) * np.exp(1j * residual_phase)
    denom = Cxy_inv*Cxy_inv - 4*Cxx_inv*Cyy_inv
    Cxx = Cyy_inv/denom
    Cyy = Cxx_inv/denom
    Cxy = -Cxy_inv/denom
    phase_quad = Cxx*xx*xx + Cyy*yy*yy + Cxy*xx*yy

    if Cxx_inv != 0:
        alpha_inv = 1j * np.sqrt(-1/denom + 0j)
    else:
        alpha_inv = 1/Cxy_inv

    # compute inverse semi analytical FFT
    fft_c = scipy.fft.ifft2(U)
    c = scipy.fft.ifftshift(fft_c)
    U = c * np.exp(1j * phase_quad)
    E = alpha_inv * scipy.fft.fft2(scipy.fft.fftshift(U)) * np.exp(1j * phase_quad)

    return E

def findFWHM(powerField, dx):
    maxVal = np.max(powerField)
    half_max = maxVal / 2.0
    maxIndex = np.array(np.unravel_index(np.argmax(powerField), powerField.shape))
    # maxIndex = (512, 512)
    left_idx = np.array(np.unravel_index(np.argmax(powerField > half_max), powerField.shape))
    fwhm = (np.sqrt(np.sum((maxIndex - left_idx)**2)))*2 * dx
    return fwhm, maxVal