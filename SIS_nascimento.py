import numpy
import matplotlib.pyplot as plt


def SIS_nascimento(N, I0, dt, beta, gamma, mu, nmax):
    S = numpy.zeros(nmax)
    I = numpy.zeros(nmax)
    S[0] = N - I0
    I[0] = I0
    for n in range(1, nmax):
        S[n] = S[n-1]*(1 - dt*beta/N*I[n-1]) + dt * \
            gamma*I[n-1] + mu*dt*(N - S[n-1])
        I[n] = I[n-1]*(1 + dt*(beta/N*S[n-1])) - dt*gamma*I[n-1] - dt*mu*I[n-1]
    return S, I


def grafico_SIS_nascimento():
    N = 1000                   # Tamanho população
    nmax = 1000                # quantidade de amostras
    I0 = 1                     # qtd. infectados  p/ n = 0
    beta = 1/7                 # taxa de transmissão / dia
    gamma = 0.03               # taxa de recuperação / dia
    mu = 0.02                  # taxa de natalidade e mortalidade
    dt = 0.1/beta              # período de amostragem
    S, I = SIS_nascimento(N, I0, dt, beta, gamma, mu, nmax)
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


grafico_SIS_nascimento()
