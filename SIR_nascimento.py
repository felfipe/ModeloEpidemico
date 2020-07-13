import numpy
import matplotlib.pyplot as plt


def SIR_nascimento(S0, I0, R0, T, beta, gamma, mu, nmax):
    S = numpy.zeros(nmax)
    I = numpy.zeros(nmax)
    R = numpy.zeros(nmax)
    S[0] = S0
    I[0] = I0
    R[0] = R0
    N = S0 + I0 + R0
    for n in range(1, nmax):
        S[n] = S[n-1]*(1 - T*beta/N*I[n-1]) + T*mu*(N - S[n-1])
        I[n] = I[n-1]*(1 + T*(beta/N*S[n-1] - gamma - mu))
        R[n] = R[n-1] + T*gamma*I[n-1] - T*mu*R[n-1]
    return S, I, R


def figura():
    N = 1000                   # Tamanho população
    nmax = 1000                # quantidade de amostras
    I0 = 1                    # qtd. infectados  p/ n = 0
    R0 = 0                     # qtd. removidos   p/ n = 0
    S0 = N - I0 - R0           # qtd. suscetíveis p/ n = 0
    beta = 1/7                 # taxa de transmissão / dia
    gamma = 0.03               # taxa de recuperação / dia
    mu = 0.01                  # taxa de natalidade e mortalidade
    dt = 0.1/beta              # período de amostragem
    Ttotal = numpy.linspace(0, dt*nmax, nmax)
    S, I, R = SIR_nascimento(S0, I0, R0, dt, beta, gamma, mu, nmax)
    Ttotal = numpy.linspace(0, dt*nmax, nmax)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(Ttotal, S, 'g--', label='Suscetíveis')
    ax.plot(Ttotal, I, 'r--', label='Infectados')
    ax.plot(Ttotal, R, 'b--', label='Removidos')
    ax.set_ylabel('Pessoas infectadas')
    ax.set_xlabel('Tempo (dias)')
    ax.grid()
    ax.legend()
    plt.xlim(0, 400)
    plt.show()


figura()
