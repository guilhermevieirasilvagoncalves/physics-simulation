from flask import Flask, render_template, request, redirect, session, flash
import numpy as np
from sympy import *
import sympy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import base64
from io import BytesIO 
from matplotlib.figure import Figure
import shutil
app = Flask(__name__) #import name

@app.route('/')
def index():
    return render_template('index.html')
    
#----------------- RENDER TELA 1 -----------------

@app.route('/tela1', methods=['GET', 'POST'])
def tela1():
    if request.method == 'GET':
        return render_template('tela1.html')
    else:
        retorno = []
        L = np.double(request.form['L'])
        ni = int(request.form['ni'])
        nf = int(request.form['nf'])
        a = np.double(request.form['a'])
        b = np.double(request.form['b'])
        Lf = L

        A = np.sqrt(2 / L)
        k = nf*np.pi / L
        retornoNF = ("nf = {}, ψ(x) = {} * sen({})·x".format(nf,np.format_float_scientific(A, precision = 2, exp_digits = 1),np.format_float_scientific(k, precision = 2, exp_digits = 1)))
        k = ni*np.pi / L
        retornoNI = ("ni = {}, ψ(x) = {} * sen({})·x".format(ni,np.format_float_scientific(A, precision = 2, exp_digits = 1),np.format_float_scientific(k, precision = 2, exp_digits = 1)))

        # plot dos gráficos das funções de onda para o nível inicial e final da partícula
        def gx(x):
            return np.sqrt(2/L)*np.sin(x*nf*np.pi/L)
        def g(x):
            return np.sqrt(2/L)*np.sin(x*ni*np.pi/L)
        #fig = Figure()
        def grafico1():
            x = np.linspace(0,L,50)
            plot1 = plt.subplot2grid((2, 2), (0, 0), colspan=2) 
            plot2 = plt.subplot2grid((2, 2), (1, 0), rowspan=2, colspan=2) 
            plot1.plot(x, gx(x)) 
            plot1.set_title('nf') 
            plot2.plot(x, g(x)) 
            plot2.set_title('ni')
            #buf = BytesIO()
            #plt.savefig(buf, format="png")
            plt.tight_layout()
            fig = plt.gcf() 
            fig.savefig('static/grafico1.png', format='png')
    
        #def plotgraf1():
            #return f"<img src='/static/images/new_plot.png'/>"    
            #data = base64.b64encode(buf.getbuffer()).decode("ascii")
            #return render_template('untitled1.html', name = 'new_plot', url ='/static/images/new_plot.png')
        
        Enf = pow(nf,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))
        Eni = pow(ni,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))
        retornoEnfp = ("E{} = {}J ou {}eV".format(nf,np.format_float_scientific(Enf, precision = 2, exp_digits = 1),np.format_float_scientific(Enf/1.602E-19, precision = 4, exp_digits = 1)))
        retornoEnip = ("E{} = {}J ou {}eV".format(ni,np.format_float_scientific(Eni, precision = 2, exp_digits = 1),np.format_float_scientific(Eni/1.602E-19, precision = 4, exp_digits = 1)))
        v = np.sqrt(2*Enf/1.67E-27)
        retornoNfVp = ("v{} = {:.2f}m/s".format(nf,v))
        v = np.sqrt(2*Eni/1.67E-27)
        retornoNiVp = ("v{} = {:.2f}m/s".format(ni,v))

        Efoton = np.subtract(Enf/1.602E-19,Eni/1.602E-19)
        retornoEfoton = ("Efóton = {:.4f}eV".format(Efoton))

        Enf = pow(nf,2)*pow(6.626E-34,2)/(8*9.11E-31*pow(L,2))
        Eni = pow(ni,2)*pow(6.626E-34,2)/(8*9.11E-31*pow(L,2))
        v = np.sqrt(2*Enf/9.11E-31)
        retornoNfVe = ("v{} = {:.2f}m/s".format(nf,v))
        v = np.sqrt(2*Eni/9.11E-31)
        retornoNiVe = ("v{} = {:.2f}m/s".format(ni,v))

        retornoEnfe = ("E{} = {}J ou {}eV".format(nf,np.format_float_scientific(Enf, precision = 2, exp_digits = 1),np.format_float_scientific(Enf/1.602E-19, precision = 4, exp_digits = 1)))
        retornoEnie = ("E{} = {}J ou {}eV".format(ni,np.format_float_scientific(Eni, precision = 2, exp_digits = 1),np.format_float_scientific(Eni/1.602E-19, precision = 4, exp_digits = 1)))

        CompFoton = 4.136E-15 * 3E8/ Efoton
        retornoCompFoton = ("λ = {}m".format(np.format_float_scientific(CompFoton, precision = 2, exp_digits = 1)))

        ldebroglie = 2*L/nf
        retornobroglieNf = ("λ De Broglie de {} = {}m".format(nf,np.format_float_scientific(ldebroglie, precision = 2, exp_digits = 1)))
        ldebroglie = 2*L/ni
        retornobrodlieNi = ("λ De Broglie de {} = {}m".format(ni,np.format_float_scientific(ldebroglie, precision = 2, exp_digits = 1)))


        L = round(L/1E-9, 3)
        thetai = nf*sp.pi*a/L
        thetaf = nf*sp.pi*b/L
        x = Symbol("x")
        p = sp.integrate((2/nf*sp.pi)*(pow(sp.sin(x),2)), (x, thetai, thetaf))
        retornoPNf =("n = {}, P({}≤x≤{}) = {}%".format(nf,thetai,thetaf,np.format_float_scientific(p*10, precision = 2, exp_digits = 1)))
        thetai = ni*sp.pi*a/L
        thetaf = ni*sp.pi*b/L
        p = sp.integrate((2/ni*sp.pi)*(pow(sp.sin(x),2)), (x, thetai, thetaf))
        RetornoPNi = ("n = {}, P({}≤x≤{}) = {}%".format(ni,thetai,thetaf,np.format_float_scientific(p*10, precision = 2, exp_digits = 1)))

        retorno.append(retornoNF)
        retorno.append(retornoNI)
        retorno.append(retornoEnfp)
        retorno.append(retornoEnip)
        retorno.append(retornoEnfe)
        retorno.append(retornoEnie)
        retorno.append(retornoEfoton)
        retorno.append(retornoCompFoton)
        retorno.append(retornoNfVp)
        retorno.append(retornoNiVp)
        retorno.append(retornoNfVe)
        retorno.append(retornoNiVe)
        retorno.append(retornobroglieNf)
        retorno.append(retornobrodlieNi)
        retorno.append(RetornoPNi)
        retorno.append(retornoPNf)
        #retorno.append(plotgraf1())

        #plot dos gráficos da função distribuição de probabilidade para o nível inicial e final da partícula
        x = np.linspace(0,L,50)
        def fx(x):
            return (2/L)*pow(np.sin(x*ni*np.pi/L),2)

        def f(x):
            return (2/L)*pow(np.sin(x*nf*np.pi/L),2)
        
        def grafico2():
            plot3 = plt.subplot2grid((2, 2), (0, 0), colspan=2) 
            plot4 = plt.subplot2grid((2, 2), (1, 0), rowspan=2, colspan=2) 

            plot3.plot(x, fx(x)) 
            plot3.set_title('ni') 
            plot4.plot(x, f(x)) 
            plot4.set_title('nf') 
            plt.tight_layout() 
            #fig = plt.gcf() 
            plt.savefig('static/grafico2.png', format='png')

        def grafico3(L):
            print(L)
            E1 = round(pow(1,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))/1.602E-19,3)
            E2 = round(pow(2,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))/1.602E-19,3)
            E3 = round(pow(3,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))/1.602E-19,3)
            E4 = round(pow(4,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))/1.602E-19,3)
            E5 = round(pow(5,2)*pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))/1.602E-19,3)
            lin1 = [E1,E1,E1,E1,E1]
            lin2 = [E2,E2,E2,E2,E2]
            lin3 = [E3,E3,E3,E3,E3]
            lin4 = [E4,E4,E4,E4,E4]
            lin5 = [E5,E5,E5,E5,E5]

            sub1 = "E1 = {}eV".format(E1)
            sub2 = "E2 = {}eV".format(E2)
            sub3 = "E3 = {}eV".format(E3)
            sub4 = "E4 = {}eV".format(E4)
            sub5 = "E5 = {}eV".format(E5)

            plt.plot(lin1, label = sub1)
            plt.plot(lin2, label = sub2)
            plt.plot(lin3, label = sub3)
            plt.plot(lin4, label = sub4)
            plt.plot(lin5, label = sub5)
            plt.legend() 
            #fig = plt.gcf() 
            plt.savefig('static/grafico3.png', format='png')
        # plotgraf1()
        # plotgraf2()
        # plotgraf3(Lf)
        grafico3(Lf)
        grafico1()
        grafico2()
        return render_template('tela1-output.html', retorno= retorno)

#----------------- RENDER TELA 2 -----------------

@app.route('/tela2', methods=['GET', 'POST'])
def tela2():
    if request.method == 'GET':
        return render_template('tela2.html')
    else:
        retorno = []
        A = np.double(request.form['A'])
        k = np.double(request.form['k'])
        #P = (request.form['P'])
        #E = (request.form['E'])
    
        L = 2/pow(A,2);
        n = (k*L)/np.pi;
        Efund = pow(6.626E-34,2)/(8*9.11E-31*pow(L,2))
        v = np.sqrt(2*Efund/9.11E-31)

        resL = ("L = " + np.format_float_scientific(L, precision = 2, exp_digits = 1) + "m or " + np.format_float_scientific(L*1000000000, precision = 2, exp_digits = 1) + "nm")
        resN = ("n = " + np.format_float_scientific(n, precision = 2, exp_digits = 1))
        resE1 = ("E1 = " + np.format_float_scientific(Efund, precision = 2, exp_digits = 1) + "J or "+ np.format_float_scientific(Efund/1.602E-19, precision = 2, exp_digits = 1) + "eV")
        resV = ("v = " + np.format_float_scientific(v, precision = 2, exp_digits = 1) + "m/s")
        retorno.append(resL)
        retorno.append(resN)
        retorno.append(resE1)
        retorno.append(resV)

        L = 2/pow(A,2);
        n = (k*L)/np.pi;
        Efund = pow(6.626E-34,2)/(8*1.67E-27*pow(L,2))
        v = np.sqrt(2*Efund/1.67E-27)

        resL = ("L = " + np.format_float_scientific(L, precision = 2, exp_digits = 1) + "m or " + np.format_float_scientific(L*1000000000, precision = 2, exp_digits = 1) + "nm")
        resN = ("n = " + np.format_float_scientific(n, precision = 2, exp_digits = 1))
        resE1 = ("E1 = " + np.format_float_scientific(Efund, precision = 2, exp_digits = 1) + "J or "+ np.format_float_scientific(Efund/1.602E-19, precision = 2, exp_digits = 1) + "eV")
        resV = ("v = " + np.format_float_scientific(v, precision = 2, exp_digits = 1) + "m/s")
        retorno.append(resL)
        retorno.append(resN)
        retorno.append(resE1)
        retorno.append(resV)

        return render_template('tela2-output.html', retorno=retorno)

#----------------- RENDER TELA 3 -----------------

@app.route('/tela3', methods=['GET', 'POST'])
def tela3():
    if request.method == 'GET':
        return render_template('tela3.html')
    else:
        retorno = []
        L = np.double(request.form['L'])
        nxi = int(request.form['nxi'])
        nyi = int(request.form['nyi'])
        nxf = int(request.form['nxf'])
        nyf= int(request.form['nyf'])

        Efund = ((1**2 + 1**2)*(pow(6.626E-34,2))/(8*9.11E-3*pow(L,2)))/1.602E-19
        Ei = ((nxi**2 + nyi**2)*(pow(6.626E-34,2))/(8*9.11E-31*pow(L,2)))
        Ef = ((nxf**2 + nyf**2)*(pow(6.626E-34,2))/(8*9.11E-31*pow(L,2)))
        li = 6.626E-34/np.sqrt(2*9.11E-31*Ei)
        lf = 6.626E-34/np.sqrt(2*9.11E-31*Ef)

        resEFund = ("{}eV".format(np.format_float_scientific(Efund, precision = 2, exp_digits = 1)))
        resEI = ("{}eV".format(np.format_float_scientific(Ei/1.602E-19, precision = 2, exp_digits = 1)))
        resEf = ("{}eV".format(np.format_float_scientific(Ef/1.602E-19, precision = 2, exp_digits = 1)))
        resLI = ("{}m".format(np.format_float_scientific(li, precision = 2, exp_digits = 1)))
        resLF = ("{}m".format(np.format_float_scientific(lf, precision = 2, exp_digits = 1)))
        retorno.append(resEFund)
        retorno.append(resEI)
        retorno.append(resEf)
        retorno.append(resLI)
        retorno.append(resLF)

        Efund = ((1**2 + 1**2)*(pow(6.626E-34,2))/(8*1.67E-27*pow(L,2)))/1.602E-19
        Ei = ((nxi**2 + nyi**2)*(pow(6.626E-34,2))/(8*1.67E-27*pow(L,2)))
        Ef = ((nxf**2 + nyf**2)*(pow(6.626E-34,2))/(8*1.67E-27*pow(L,2)))
        li = 6.626E-34/np.sqrt(2*9.11E-31*Ei)
        lf = 6.626E-34/np.sqrt(2*9.11E-31*Ef)

        resEFund = ("{}eV".format(np.format_float_scientific(Efund, precision = 2, exp_digits = 1)))
        resEI = ("{}eV".format(np.format_float_scientific(Ei/1.602E-19, precision = 2, exp_digits = 1)))
        resEf = ("{}eV".format(np.format_float_scientific(Ef/1.602E-19, precision = 2, exp_digits = 1)))
        resLI = ("{}m".format(np.format_float_scientific(li, precision = 2, exp_digits = 1)))
        resLF = ("{}m".format(np.format_float_scientific(lf, precision = 2, exp_digits = 1)))
        retorno.append(resEFund)
        retorno.append(resEI)
        retorno.append(resEf)
        retorno.append(resLI)
        retorno.append(resLF)


        return render_template('tela3-output.html', retorno= retorno) 

#----------------- RENDER TELA 4 -----------------
   
@app.route('/tela4', methods=['GET', 'POST'])
def tela4():
    if request.method == 'GET':
        return render_template('tela4.html')
    else:
        retorno = []
        L = np.double(request.form['L'])
        nxi = int(request.form['nxi'])
        nyi = int(request.form['nyi'])
        nzi = int(request.form['nzi'])
        nxf = int(request.form['nxf'])
        nyf= int(request.form['nyf'])
        nzf = int(request.form['nzf'])

        Efund = ((1**2 + 1**2 + 1**2)*(pow(6.626E-34,2))/(8*9.11E-31*pow(L,2)))/1.602E-19
        Ei = ((nxi**2 + nyi**2 + nzi**2)*(pow(6.626E-34,2))/(8*9.11E-31*pow(L,2)))
        Ef = ((nxf**2 + nyf**2 + nzf**2)*(pow(6.626E-34,2))/(8*9.11E-31*pow(L,2)))
        Efoton = Ei/1.602E-19 - Ef/1.602E-19
        lfoton = 4.136E-15*3E8/Efoton

        resEfund = ("{}eV".format(np.format_float_scientific(Efund, precision = 2, exp_digits = 1)))
        resEI = ("{}eV".format(np.format_float_scientific(Ei/1.602E-19, precision = 2, exp_digits = 1)))
        resEF =("{}eV".format(np.format_float_scientific(Ef/1.602E-19, precision = 2, exp_digits = 1)))
        resEfoton = ("{:.3f}eV".format(Efoton))
        resLfotoon = ("{}m".format(lfoton))

        retorno.append(resEfund)
        retorno.append(resEI)
        retorno.append(resEF)
        retorno.append(resEfoton)
        retorno.append(resLfotoon)

        Efund = ((1**2 + 1**2 + 1**2)*(pow(6.626E-34,2))/(8*1.67E-27*pow(L,2)))/1.602E-19
        Ei = ((nxi**2 + nyi**2 + nzi**2)*(pow(6.626E-34,2))/(8*1.67E-27*pow(L,2)))
        Ef = ((nxf**2 + nyf**2 + nzf**2)*(pow(6.626E-34,2))/(8*1.67E-27*pow(L,2)))
        Efoton = Ei/1.602E-19 - Ef/1.602E-19
        lfoton = 4.136E-15*3E8/Efoton

        resEfund = ("{}eV".format(np.format_float_scientific(Efund, precision = 2, exp_digits = 1)))
        resEI = ("{}eV".format(np.format_float_scientific(Ei/1.602E-19, precision = 2, exp_digits = 1)))
        resEF =("{}eV".format(np.format_float_scientific(Ef/1.602E-19, precision = 2, exp_digits = 1)))
        resEfoton = ("{:.3f}eV".format(Efoton))
        resLfotoon = ("{}m".format(lfoton))
        retorno.append(resEfund)
        retorno.append(resEI)
        retorno.append(resEF)
        retorno.append(resEfoton)
        retorno.append(resLfotoon)

        return render_template('tela4-output.html', retorno= retorno) 

if __name__ == "__main__":
    app.run(debug = True)
