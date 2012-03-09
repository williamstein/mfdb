
from flask import Flask, request, redirect, url_for, g, render_template

app = Flask(__name__)



@app.route("/")
def main_page():
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
