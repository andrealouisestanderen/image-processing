from flask import Flask, request, render_template
import os
import matplotlib
matplotlib.use('Agg')
import analysis

print("BAR CHART: http://127.0.0.1:5000/bar-chart?file=globalterrorism.csv&column=iyear&sort=0")
print("LINE CHART: http://127.0.0.1:5000/line-chart?file=globalterrorism.csv&year=iyear&group=region_txt")

app = Flask(__name__)

@app.route("/")
def hello():
    return render_template('home.html')

@app.route("/upload")
def uploadImage():
    