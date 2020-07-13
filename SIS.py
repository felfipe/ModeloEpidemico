import numpy
import matplotlib.pyplot as plt


def SIS(N, I0, dt, beta, gama, nmax):
    S = numpy.zeros(nmax)
    I = numpy.zeros(nmax)
    S[0] = N - I0
    I[0] = I0
    for n in range(1, nmax):
        S[n] = S[n-1]*(1 - dt*beta/N*I[n-1]) + dt*gama*I[n-1]
        I[n] = I[n-1]*(1 + dt*(beta/N*S[n-1])) - dt*gama*I[n-1]
    return S, I


def grafico_SIS():
    beta = 0.1      # taxa de transmissão
    gama = 0.01     # taxa de recuperação
    N = 1000
    dt = 0.1/beta
    nmax = 1000
    S, I = SIS(N, 1, dt, beta, gama, nmax)
    Ttotal = numpy.linspace(0, dt*nmax, nmax)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(Ttotal, S, 'g--', label='suscetíveis')
    ax.plot(Ttotal, I, 'r--', label='infectados')
    ax.set_ylabel('Pessoas infectadas')
    ax.set_xlabel('Tempo (dias)')
    ax.grid()
    ax.legend()
    plt.xlim(0, 200)
    plt.show()


grafico_SIS()
