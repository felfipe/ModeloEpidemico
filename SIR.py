import numpy
import matplotlib.pyplot as plt


def SIR(S0, I0, R0, T, beta, gamma, nmax):
    S = numpy.zeros(nmax)
    I = numpy.zeros(nmax)
    R = numpy.zeros(nmax)
    S[0] = S0
    I[0] = I0
    R[0] = R0
    N = S0 + I0 + R0
    for n in range(1, nmax):
        S[n] = S[n-1]*(1 - T*beta/N*I[n-1])
        I[n] = I[n-1]*(1 + T*(beta/N*S[n-1] - gamma))
        R[n] = R[n-1] + T*gamma*I[n-1]
    return S, I, R


def curva_infectados_SIR():
    N = 1000
    I0 = 50
    R0 = 0
    S0 = N - I0 - R0
    beta = 1/7                 # taxa de transmissão / dia
    gamma = 0.03               # taxa de recuperação / dia
    nmax = 1000
    dt = 0.1/beta              # período de amostragem
    Ttotal = numpy.linspace(0, dt*nmax, nmax)
    S, I1, R = SIR(S0, I0, R0, dt, 0.2, gamma, nmax)
    S, I2, R = SIR(S0, I0, R0, dt, 0.4, gamma, nmax)
    S, I3, R = SIR(S0, I0, R0, dt, 0.6, gamma, nmax)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(Ttotal, I1, 'g--', label='beta = 0.2')
    ax.plot(Ttotal, I2, 'r--', label='beta = 0.4')
    ax.plot(Ttotal, I3, 'b--', label='beta = 0.6')
    ax.set_ylabel('Pessoas infectadas')
    ax.set_xlabel('Tempo')
    ax.grid()
    ax.legend()
    plt.xlim(0, 150)
    plt.show()


curva_infectados_SIR()
