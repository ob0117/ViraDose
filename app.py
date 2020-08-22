from flask import Flask
from flask import render_template
from flask import request, redirect

import ViraDose

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def index():

  if request.method == "POST":

    data = request.form
    ViraDose.execute(data)
    
  return render_template("index.html")

if __name__ == "__main__":
  app.run()