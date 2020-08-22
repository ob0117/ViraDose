## ViraDose Science Fair Project

ViraDose is a project created for the 2019 Waterloo-Regional and Canada-Wide Science Fairs. Originally, it was made using only Python.
This updated version is a web app which displays the full project and makes it easier to use. 

Using the user's input data, ViraDose calculates an ideal gene therapy clinical trial starting dose. It can then display the predicted 
innate immune response, risk of toxicity, or estimated effectiveness of the dose. The idea behind this app was to make gene therapy
more viable in mainstream medicine.

Running app.py opens the main page which contains introductory information, the Canada-Wide Science Fair display board, and the actual form for 
filling out data to use the app. Using the app will open a new tab displaying the calculated dose and corresponding figures. The Python project 
uses SciPy's odeint method for solving systems of differential equations. Graphs are created with matplotlib using mpld3 to display them in the browser.
