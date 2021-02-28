# IMPORTS ----------------------------------------------------------------------------------------------------------------- #

import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.special
import math

import mpld3
from mpld3 import plugins

import webbrowser
import os


def execute(data):   # function to be called within app.py (obtains html form data)

    # SETTING VARIABLES -------------------------------------------------------------------------------------------------------- #
        
    # Variables obtained from HTML form

    animal = data["animal"]
    aav = data["aavVectors"]

    doseVal = float(data["hnstdVal"])
    doseExp = int(data["hnstdExp"])

    weight1 = float(data["scalingWeight"])
    weight2 = float(data["patientWeight"])
    safetyFactor = float(data["safetyFactor"])

    mathModel = data["mathModel"]

    genomeweightVal = float(data["genomeWeightVal"])
    genomeweightExp = int(data["genomeWeightExp"])
    genomelength = float(data["genomeLength"])
    IFN = float(data["IFNinit"])
    AVP = float(data["AVPinit"])

    t = float(data["t"])
    logc = float(data["logc"])
    c1 = float(data["c1"])
    h = float(data["h"])
    h1 = float(data["hh"])

    k1 = 0
    k2 = float(data["k2"])
    k3 = float(data["k3"]) 
    d1 = float(data["d1"])
    d2 = float(data["d2"])
    d3 = float(data["d3"])
    b1 = float(data["b1"])
    b2 = float(data["b2"])
    K1 = float(data["capk1"]) 
    K2 = float(data["capk2"])
    n1 = float(data["n1"])
    n2 = float(data["n2"]) 

    # Setting all other variables for use

    animalDict = {
        "Mouse" : 12.3 ,
        "Rat" : 6.2 ,
        "Hamster" : 7.4 ,
        "Ferret" : 5.3 ,
        "Guinea Pig" : 4.6 ,
        "Rabbit" : 3.1 ,
        "Dog" : 1.8 ,
        "Rhesus Monkey" : 3.1 ,
        "Marmoset" : 6.2 ,
        "Squirrel Monkey" : 5.3 ,
        "Baboon" : 1.8 ,
        "Micro Pig" : 1.4 ,
        "Mini Pig" : 1.1
    }

    HNSTD = doseVal * (10 ** doseExp)
    animalFactor = animalDict[animal]
    genomeweight = genomeweightVal * (10 ** genomeweightExp)

    e = math.e



    # CALCULATING FLAT DOSE ---------------------------------------------------------------------------------------------------- #

    HED =  HNSTD / animalFactor # human equivalent dose

    if safetyFactor == 0:  # applying safety factor
        safe = HED
    else:
        safe = HED / safetyFactor

    Flat = safe * weight1  # scaling dose by scaling weight




    # IF DISPLAYING MATH MODEL 1 ----------------------------------------------------------------------------------------------- #

    if mathModel == "View innate immune response" :

        # a rough conversion of the dose to mol/L concentration (to fit the units in the model)

        if genomeweight == 5.1e-18 and genomelength == 4700:
            genomeweight = 5.10185 * (10 ** (-18))
            Da = 650
        else:
            Da = genomeweight / (1.67 * (10 ** (-24)))

        g = genomeweight * Flat

        gMol = Da

        mol = g / gMol

        L = weight2 * 0.07 * (0.55 - (0.55 * (0.1) * (0.07)))
        Vm = mol / L

        # Defining system of differential equations

        def immunity(y, t):

            V = y[0]
            I = y[1]
            A = y[2]

            dVdt = k1 * V * ((b1 * (K1 ** n1)) / ((K1 ** n1) + (A ** n1))) - d1 * V
            dIdt = k2 * V + ((b2 * (I ** n2)) / ((K2 ** n2) + (I ** n2))) - d2 * I
            dAdt = k3 * I - d3 * A

            return [dVdt, dIdt, dAdt]

        y0 = [Vm, IFN, AVP]   # initial values
        t = np.linspace(0, 200, 500)  # defining time points for input (range and number of points)

        # Solving the system based on initial values + setting up graph points

        y = odeint(immunity, y0, t)   # odeint numerically solves system given the equations, initial values and time points
        
        # SET UP GRAPH

        V = y[:, 0]
        I = y[:, 1]
        A = y[:, 2]

        fig = plt.figure(figsize = (10,6))

        Vpoints = plt.plot(t, V, label='Viral mRNA', marker="o", markersize="1")
        Ipoints = plt.plot(t, I, label='IFN', marker="o", markersize="1")
        Apoints = plt.plot(t, A, label='AVP', marker="o", markersize="1")
        plt.legend(loc='upper right')
        plt.xlabel('Time Post-Infection (hours)')
        plt.ylabel('Concentration (mol/L)')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.title('Model of Predicted Immune Response')

        # mpld3 tooltip plugin makes points show up on web graph
        # points2 = plt.plot(xvaluesm, toxfunction, marker="o", markersize="1", color="cornflowerblue")
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(Vpoints[0]))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(Ipoints[0]))
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(Apoints[0]))


        graphString = mpld3.fig_to_html(fig)  # turns graph into HTML code to be embedded on webpage

        # Write graph and formatting to HTML file

        htmlFile = open("graph.html","w")

        htmlFile.write("<!DOCTYPE html> <html>")
        htmlFile.write("<link rel='stylesheet' href='static/AppStyle.css'>") 
        htmlFile.write("<head>   <meta charset = 'UTF-8' />   <meta name='viewport' content='width=device-width initial-scale=1'>")
        htmlFile.write("<meta name='viewport' content='width=device-width initial-scale=1'>  <title>ViraDose Project Display</title>")
        htmlFile.write("<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css' integrity='sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm' crossorigin='anonymous'>")
        htmlFile.write("<script src='https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js'></script>")
        htmlFile.write("<script src='https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.16.0/umd/popper.min.js'></script>")
        htmlFile.write("<script src='https://maxcdn.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.min.js'></script>")

        htmlFile.write("<body style='height:auto; text-align:center'> <h1 class='graph-heading'> Recommended starting dose: " + "{:.2e}".format(Flat) +  " genome copies </h1>")

        htmlFile.write(graphString)

        htmlFile.write("<br><br><br><h5> Input values used: </h5> <br>")
        htmlFile.write("<h6> <b>HNSTD</b>: " + "{:.2e}".format(HNSTD) +  " genome copies per kilogram </h6> <br>")
        htmlFile.write("<h6> <b>AAV vector</b>: " + aav + "</h6> <br>")
        htmlFile.write("<h6> <b>Animal</b>: " + animal + "</h6> <br>")
        htmlFile.write("<h6> <b>Scaling weight</b>: " + str(weight1) + "</h6> <br>")
        htmlFile.write("<h6> <b>Safety factor</b>: " + str(safetyFactor) + "</h6> <br>")
        htmlFile.write("<h6> <b>Patient's weight</b>: " + str(weight2) + "</h6> <br>")
        htmlFile.write("<h6> <b>Genome weight</b>: " + "{:.2e}".format(genomeweight) + "</h6> <br>")
        htmlFile.write("<h6> <b>Genome length</b>: " + str(genomelength) + "</h6> <br>")
        htmlFile.write("<h6> <b>Initial IFN</b>: " + str(IFN) + "</h6> <br>")
        htmlFile.write("<h6> <b>Initial AVP</b>: " + str(AVP) + "</h6> <br>")
        htmlFile.write("<h6> <b>k2</b>: " + str(k2) + "</h6> <br>")
        htmlFile.write("<h6> <b>k3</b>: " + str(k3) + "</h6> <br>")
        htmlFile.write("<h6> <b>d1</b>: " + str(d1) + "</h6> <br>")
        htmlFile.write("<h6> <b>d2</b>: " + str(d2) + "</h6> <br>")
        htmlFile.write("<h6> <b>d3</b>: " + str(d3) + "</h6> <br>")
        htmlFile.write("<h6> <b>b1</b>: " + str(b1) + "</h6> <br>")
        htmlFile.write("<h6> <b>b2</b>: " + str(b2) + "</h6> <br>")
        htmlFile.write("<h6> <b>K1</b>: " + str(K1) + "</h6> <br>")
        htmlFile.write("<h6> <b>K2</b>: " + str(K2) + "</h6> <br>")
        htmlFile.write("<h6> <b>n1</b>: " + str(n1) + "</h6> <br>")
        htmlFile.write("<h6> <b>n2</b>: " + str(n2) + "</h6> <br><br><br><br>")
        htmlFile.write("</body></html>")

        htmlFile.close()

        webbrowser.open("graph.html")  # open HTML file




    # IF DISPLAYING MATH MODEL 2 --------------------------------------------------------------------------------------------------------- #

    else :

        livercells = (139e6) * 1500  # number of hepatocytes per gram * grams of liver
        # ^this represents target cells, since liver is common target organ (for now)

        m = Flat / livercells  # input dose per target cell (is what this model uses)
        
        # Function for gene expression (takes number of cells as input)
        def expression(k):

            c = 10 ** logc

            v = t / (1 + ((c / m) ** h))

            P = ((v ** k) * (e ** (-v))) / (scipy.special.factorial(k))

            return P

        # Function for toxicity (takes dose as input)
        def toxicity(m):

            Pd = 1 / (1 + ((c1 / m) ** h1))

            return Pd
        
        fig = plt.figure(figsize = (10,6))

        # PLOTTING EXPRESSION

        xvaluesk = np.linspace(0, 10, 1000)
        expfunction = expression(xvaluesk)

        #plt.figure(1)

        # mpld3 tooltip plugin makes points show up on web graph
        points = plt.plot(xvaluesk, expfunction, marker="o", markersize="1")
        mpld3.plugins.connect(fig, mpld3.plugins.PointLabelTooltip(points[0]))

        plt.xlabel('Number of Gene-Expressing Vectors in Cell')
        plt.ylabel('Probability of Random Cell Containing That Number')

        plt.title('Predicted Transgene Expression')

        # PLOTTING TOXICITY

        fig2 = plt.figure(figsize = (10,6))

        xvaluesm = np.linspace(0, (m * 1000), 1000)
        toxfunction = toxicity(xvaluesm)

        #plt.figure(2)

        plt.plot(xvaluesm, toxfunction, label="Theoretical Dose-Dependent Toxicity")

        plt.xlabel("Input Dose per Target Cell")
        plt.ylabel("Proportion of Apoptotic Cells")
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

        # Adding points to show where our current dose is on the input axis, as well as where ten-fold increases lie
        plt.plot(m, toxicity(m), marker='o', markersize=6, color="black", label="Current Input Dose")
        plt.plot(m * 10, toxicity(m * 10), marker='o', markersize=6, color="red", label="Ten-fold Dose Increases")
        plt.plot(m * 100, toxicity(m * 100), marker='o', markersize=6, color="red")
        plt.plot(m * 1000, toxicity(m * 1000), marker='o', markersize=6, color="red")

        # mpld3 tooltip plugin makes points show up on web graph
        points2 = plt.plot(xvaluesm, toxfunction, marker="o", markersize="1", color="cornflowerblue")
        mpld3.plugins.connect(fig2, mpld3.plugins.PointLabelTooltip(points2[0]))

        plt.legend(loc="upper right")

        plt.title("Predicted Cytotoxicity")

        graphString1 = mpld3.fig_to_html(fig2)  # turning graphs into HTML code to be embedded on webpage
        graphString2 = mpld3.fig_to_html(fig)

        # Writing HTML file

        htmlFile = open("graph.html","w")

        htmlFile.write("<!DOCTYPE html> <html>")
        htmlFile.write("<link rel='stylesheet' href='static/AppStyle.css'>") 
        htmlFile.write("<head>   <meta charset = 'UTF-8' />   <meta name='viewport' content='width=device-width initial-scale=1'>")
        htmlFile.write("<meta name='viewport' content='width=device-width initial-scale=1'>  <title>ViraDose Project Display</title>")
        htmlFile.write("<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css' integrity='sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm' crossorigin='anonymous'>")
        htmlFile.write("<script src='https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js'></script>")
        htmlFile.write("<script src='https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.16.0/umd/popper.min.js'></script>")
        htmlFile.write("<script src='https://maxcdn.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.min.js'></script>")

        htmlFile.write("<body style='height:auto; text-align:center'> <h1 class='graph-heading'> Recommended starting dose: " + "{:.2e}".format(Flat) +  " genome copies </h1>")

        htmlFile.write(graphString1)
        htmlFile.write(graphString2)

        htmlFile.write("<br><br><br><h5> Input values used: </h5> <br>")
        htmlFile.write("<h6> <b>HNSTD</b>: " + "{:.2e}".format(HNSTD) +  " genome copies per kilogram </h6> <br>")
        htmlFile.write("<h6> <b>AAV vector</b>: " + aav + "</h6> <br>")
        htmlFile.write("<h6> <b>Animal</b>: " + animal + "</h6> <br>")
        htmlFile.write("<h6> <b>Scaling weight</b>: " + str(weight1) + "</h6> <br>")
        htmlFile.write("<h6> <b>Safety factor</b>: " + str(safetyFactor) + "</h6> <br>")
        htmlFile.write("<h6> <b>Patient's weight</b>: " + str(weight2) + "</h6> <br>")
        htmlFile.write("<h6> <b>Upper limit</b>: " + str(t) + "</h6> <br>")
        htmlFile.write("<h6> <b>log EC50</b>: " + str(logc) + "</h6> <br>")
        htmlFile.write("<h6> <b>Half-maximal toxicity</b>: " + str(c1) + "</h6> <br>")
        htmlFile.write("<h6> <b>Hill Coefficient (gene expression)</b>: " + str(h) + "</h6> <br>")
        htmlFile.write("<h6> <b>Hill Coefficient (toxicity)</b>: " + str(h1) + "</h6> <br><br><br><br>")

        htmlFile.write("</body> </html>")

        htmlFile.close()

        webbrowser.open("graph.html")  # opening HTML file