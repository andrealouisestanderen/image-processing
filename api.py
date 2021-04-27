from flask import Flask, render_template, flash, request, redirect, url_for, send_from_directory
from werkzeug.utils import secure_filename
from datetime import datetime
import os
import subprocess
import time

UPLOAD_FOLDER = './Static'
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@app.route('/', methods=['GET', 'POST'])
def upload_file():
    # Delete older processed images
    for filename in os.listdir('static/'):
        if filename.startswith('blurry'):  # not to remove other images
            os.remove('static/' + filename)
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit a empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

            # can add if-else for selecting filter type
            return redirect(url_for('processed_file',
                                    filename=filename))
    return render_template('home.html')


@app.route('/processed/<filename>')
def processed_file(filename):

    newFilename = "blurry" + str(time.time()) + ".png"
    subprocess.call(['./processImage.sh', filename, newFilename])

    return render_template('processed.html', file=newFilename)


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
