# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 12:55:43 2017

@author: pbarros
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from graph import graph

#Locale settings (para gráficos com vírgulas)
#import locale
# Set to German locale to get comma decimal separater

#locale.setlocale(locale.LC_NUMERIC, "de_DE")
#plt.rcdefaults()
# Tell matplotlib to use the locale we set above
#plt.rcParams['axes.formatter.use_locale'] = True



V = 31
wmin = 0.01
wmax = 10
dw = 0.01

w = np.arange(wmin,wmax+dw,dw)
T = 2*np.pi/w
f = 1/T
C = 0.002

Sd = 916700/(2*np.pi)*C*w/(1+(191*w/V)**2)**(4/3)
Sh = 1146*0.002*V/(2+(286*w/V)**2)**(5/6)

So = np.zeros(len(w))
x = 1.592*w/V
for i in range(len(w)):
    if w[i]/V < 0.001885:
        Fg = 583*x[i]
    elif w[i]/V < 0.0628:
        Fg = 420*x[i]**0.7/(1+x[i]**0.35)**11.5
    else:
        Fg = 838*x[i]/(1+x[i]**0.35)**11.5
    So[i] = (750+69*V)*10**-6*V**2*Fg/w[i]
    
 

fig = plt.figure(figsize=(8,4),facecolor='white')    
ax = fig.add_subplot(111)
graph(fig, ax, w, So, c = 'r', n = 'Ochi-Shin', semilogx = True, xlim = [0.01,10], title = 'Espectros de Vento', lab_x = 'Frequência (rad/s)', lab_y = 'Densidade Espectral (m²/s)')
graph(fig, ax, w, Sh, c = 'b', n = 'Harris', semilogx = True)          
graph(fig, ax, w, Sd, c = 'k', l = '--', n = "Davenport", semilogx = True)




Fl = 1e-10
Ft = 2381800
Mz = -125759000
F1 = Ft



#Espectro de Forças Harris
m0h = sum(Sh*dw)
F2h = F1/V**2*m0h
Fh = F1 + F2h
S = Sh
S1 = (2*F1/V)**2*S
S2 = np.zeros(len(S))
for i in range(len(S)):
    Sj = np.zeros(len(S))
    Sj[:len(S)-i]=S[i:]
    S2[i] = 8*(F1/V**2)**2*sum(S*Sj)*dw
Sfh = S1 + S2




#Gráfico de Espectro de Forças Harris
fig2 = plt.figure(figsize=(8,4),facecolor='white')    
ax2 = fig2.add_subplot(111)
graph(fig2, ax2, w, Sfh, c = 'b', n = 'Harris', semilogx = True, xlim = [0.01,10], title = 'Espectros de Forças', lab_x = 'Frequência (rad/s)', lab_y = 'Densidade Espectral (N²s)')
graph(fig2, ax2, w, S1, c = 'r', n = 'S1', semilogx = True)          
graph(fig2, ax2, w, S2, c = 'g', n = 'S2', semilogx = True)


#Espectro de Forças Ochi-Shin
m0o = sum(So*dw)
F2o = F1/V**2*m0o
Fo = F1 + F2o
S = So
S1o = (2*F1/V)**2*S
S2o = np.zeros(len(S))
for i in range(len(S)):
    Sj = np.zeros(len(S))
    Sj[:len(S)-i]=S[i:]
    S2o[i] = 8*(F1/V**2)**2*sum(S*Sj)*dw
Sfo = S1o + S2o


#Gráfico de Espectro de Forças Ochi-Shin
fig3 = plt.figure(figsize=(8,4),facecolor='white')    
ax3 = fig3.add_subplot(111)
graph(fig3, ax3, w, Sfo, c = 'b', n = 'Ochi-Shin', semilogx = True, xlim = [0.01,10], title = 'Espectros de Forças', lab_x = 'Frequência (rad/s)', lab_y = 'Densidade Espectral (N²s)')
graph(fig3, ax3, w, S1o, c = 'r', n = 'S1', semilogx = True)          
graph(fig3, ax3, w, S2o, c = 'g', n = 'S2', semilogx = True)


#Vetor de angulos de fase
fi = 2*np.pi*np.random.random(len(S))

#Séries temporais (Harris)
Vhw = np.sqrt(2*Sh*dw)
Fhw = np.sqrt(2*Sfh*dw)
dt1 = 0.25
tempo = np.arange(40000)*dt1
Fht = np.zeros(len(tempo))
#Fht2 = np.zeros(len(tempo))
Vht = np.zeros(len(tempo))
for i, t in enumerate(tempo):
    Vht[i] = V + sum(Vhw*np.cos(w*t+fi))     
    #Fht2[i] = Fh + sum(Fhw*np.cos(w*t+fi))
    Fht[i] = F1/V**2*Vht[i]*abs(Vht[i])
mi_fh = np.average(Fht)

"""
fig4 = plt.figure(figsize=(8,4),facecolor='white')    
ax4 = fig4.add_subplot(111)
graph(fig4, ax4, tempo, Fht/1000, ylim = [0.,1.1*max(Fht)/1000], title = 'Série Temporal da Força do Vento (Harris)', lab_x = 'Tempo (s)', lab_y = 'Força (kN)')

fig5 = plt.figure(figsize=(8,4),facecolor='white')    
ax5 = fig5.add_subplot(111)
graph(fig5, ax5, tempo, Vht, ylim = [0.,1.1*max(Vht)], title = 'Série Temporal do Vento (Harris)', lab_x = 'Tempo (s)', lab_y = 'Velocidade (m/s)')
"""

#Séries temporais (Ochi-Shin)
Vow = np.sqrt(2*So*dw)
Fow = np.sqrt(2*Sfo*dw)
Fot = np.zeros(len(tempo))
Fot2 = np.zeros(len(tempo))
Vot = np.zeros(len(tempo))
for i, t in enumerate(tempo):
    Vot[i] = V + sum(Vow*np.cos(w*t+fi))     
    Fot2[i] = Fo + sum(Fow*np.cos(w*t+fi))
    Fot[i] = F1/V**2*Vot[i]*abs(Vot[i])
mi_fo = np.average(Fot)
mi_fo2 = np.average(Fot2)
a = np.var(Fot)
b = np.var(Fot2)
print(a/b)


fig6 = plt.figure(figsize=(8,4),facecolor='white')    
ax6 = fig6.add_subplot(111)
graph(fig6, ax6, tempo, Fot/1000, xlim = [0.,500], ylim = [0.,1.1*max(Fot)/1000], title = 'Série Temporal da Força do Vento (Ochi-Shin)', lab_x = 'Tempo (s)', lab_y = 'Força (kN)')

fig7 = plt.figure(figsize=(8,4),facecolor='white')    
ax7 = fig7.add_subplot(111)
graph(fig7, ax7, tempo, Vot, xlim = [0.,500], ylim = [0.,1.1*max(Vot)], title = 'Série Temporal do Vento (Ochi-Shin)', lab_x = 'Tempo (s)', lab_y = 'Velocidade (m/s)')

fig7_2 = plt.figure(figsize=(8,4),facecolor='white')    
ax7_2 = fig7_2.add_subplot(111)
graph(fig7_2, ax7_2, tempo, Fot2/1000, xlim = [0.,500], ylim = [-.1*max(Fot)/1000,1.1*max(Fot)/1000], title = 'Série Temporal da Força do Vento (Ochi-Shin)', lab_x = 'Tempo (s)', lab_y = 'Força (kN)')

#autoc = np.correlate(Fht-mi_fh,Fht-mi_fh,'full')
#autoc = np.correlate(Vht-V,Vht-V,'full')
autoc = np.correlate(Fot-mi_fo,Fot-mi_fo,'full')
#autoc = np.correlate(Vot-V,Vot-V,'full')

autoc = autoc/len(tempo)#normalizando a autocorrelação! A normalização é a divisão pelo número de pontos considerados, assim ela não vai dependendo. Se não fizer isso, o fato da amostra ter mais pontos com dt menor leva a um valor de correlaçao maior para um mesmo intervalo de amostra.
# com isso o valor limite da autocorrelação quando dt tende a zero é uma constante para cada valor de t, i.e. o grafico não vai mudar?! pelo menos não o valor máximo...
# notar que dt/Tmax = 1/len(tempo), ou seja a integral que define 1/T*integrate(x1,x2,dt) com dt constante dá a função de correlação
tempo2 = np.arange(len(autoc))*dt1

fig8 = plt.figure(figsize=(8,4),facecolor='white')    
ax8 = fig8.add_subplot(111)
graph(fig8, ax8, tempo2-len(tempo)*dt1, autoc, title = 'Autocorrelação', lab_x = 'Tempo (s)', lab_y = 'Correlação Normalizada (N²)') 


dt2 = tempo2[1]
sp = np.fft.fft(autoc/len(tempo)) #multiplica por 1/N para normalização da FFT que não foi feita pela função do módulo FFT. Nesse caso a transformação inversa não deve ter o fator 1/N para que o ciclo seja unitário
#pegar a fórmula da FFt no manual do numpy
freq = np.fft.fftfreq(len(autoc), d=dt2)
data = pd.DataFrame({'rad/s':2*math.pi*freq,'Imaginario':sp.imag,'Real':sp.real,
                         'abs(A)*dw':np.abs(sp),'2*abs(A)':2*np.abs(sp),'phase(A)':np.angle(sp)}) 
                         #notar que as energias foram concentradas em valores discretos de frequencia quando fizemos a transformação para obter o Vw
                         #por isso, é preciso dividir por dw para obter os valores do grafico contínuo

msk = (data.loc[:,'rad/s']>0.) & (data.loc[:,'rad/s']<10)
x = np.array(data.loc[msk,'rad/s'])
y = np.array(data.loc[msk,'2*abs(A)'])

#É necessário arrumar o gráfico para plotar somente os valores do gráfico contínuo onde está presente a energia
#porque o gráfico obtido é em valores discretos, pois perdemos informaçao ao termos um dw que não tende a zero
#isso é outro indicio de que o sinal foi criado artificialmente e que perdemos informação
fig9 = plt.figure(figsize=(8,4),facecolor='white')    
ax9 = fig9.add_subplot(111)
graph(fig9, ax9, x, y, title = 'Densidade Espectral da Autocorrelação', lab_x = 'Tempo (s)', lab_y = 'Densidade Espectral (N²s)', xlim = [0.01,10], semilogx = True)

"""
fig6
fig7
fig7_2
fig8
fig9
"""