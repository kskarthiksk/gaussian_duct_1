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
    fft_c = scipy.fft.fft(E)
    c = scipy.fft.fftshift(fft_c)

    kx = 2 * np.pi * scipy.fft.fftshift(scipy.fft.fftfreq(Nx, d = dx))

    argument = k2 - kx*kx

    #Calculate the propagating and the evanescent (complex) modes
    tmp = np.sqrt(argument)
    # kz = np.where(argument >= 0, tmp, 1j*tmp)
    kz = tmp

    # propagate the angular spectrum a distance z
    E = scipy.fft.ifft(scipy.fft.ifftshift(c * np.exp(1j * kz* z)))
    # E = np.fft.ifft2(c)

    return E

def quadratic_taylor_terms(func, Nx, dx):
    Cxx = (func[Nx//2+1] - 2*func[Nx//2] + func[Nx//2-1]) / (2*dx*dx)
    return Cxx


def semi_analytical_asm(E, z, k2, Nx, dx, x):
    """
    Angular spectruum method from: https://rafael-fuente.github.io/simulating-diffraction-patterns-with-the-angular-spectrum-method-and-python.html

    Modified by Karthik Sreekumar
    """

    # ---------------------------
    # Forward semi analytical FFT
    # ---------------------------

    # computing quadratic and residual phase
    phase = np.atan2(E.imag, E.real)
    Cxx = quadratic_taylor_terms(phase, Nx, dx)
    phase_quad = Cxx*x*x
    residual_phase = phase - phase_quad

    U = np.abs(E) * np.exp(1j*residual_phase)

    # calculating analytical part of fourier transform
    kx = 2 * np.pi * scipy.fft.fftshift(scipy.fft.fftfreq(Nx, d = dx))

    Cxx_inv = -1/(4*Cxx)
    phase_quad_inv = Cxx_inv*kx*kx
    alpha = np.sqrt(1j / (2* Cxx))

    # compute forward semi analytical FFT
    fft_c = scipy.fft.fft(U)
    c = scipy.fft.fftshift(fft_c)
    U = c * np.exp(1j * phase_quad_inv / (4 * Cxx * Cxx))
    U = alpha * scipy.fft.ifft(scipy.fft.ifftshift(U)) * np.exp(1j * phase_quad_inv)

    argument = k2 - kx*kx
    E = U
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
    Cxx_inv= quadratic_taylor_terms(phase, Nx, np.abs(kx[1] - kx[0]))
    print(phase)
    phase_quad_inv = Cxx_inv*kx*kx
    residual_phase = phase - phase_quad_inv

    U = np.abs(E) * np.exp(1j * residual_phase)
    Cxx = -1/4*Cxx_inv
    phase_quad = Cxx*x*x

    alpha_inv = np.sqrt(1j / (2*Cxx_inv))

    # compute inverse semi analytical FFT
    fft_c = scipy.fft.ifft(U)
    c = scipy.fft.ifftshift(fft_c)
    U = c * np.exp(1j * phase_quad / (4*Cxx_inv * Cxx_inv))
    E = alpha_inv * scipy.fft.fft(scipy.fft.fftshift(U)) * np.exp(1j * phase_quad)

    return E

def findFWHM(powerField, dx):
    maxVal = np.max(powerField)
    half_max = maxVal / 2.0
    maxIndex = np.array(np.unravel_index(np.argmax(powerField), powerField.shape))
    # maxIndex = (512, 512)
    left_idx = np.array(np.unravel_index(np.argmax(powerField > half_max), powerField.shape))
    fwhm = (np.sqrt(np.sum((maxIndex - left_idx)**2)))*2 * dx
    return fwhm, maxVal