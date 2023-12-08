from flask import Flask, render_template, request, redirect, url_for, make_response
from validation_tool.constraints.constraints import Constraints
import analyse_input

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def validate():
    form = {'gcContentMinPercentage': 25, 
            'gcContentMaxPercentage': 65, 
            'maxHomopolymer': 5,
            'maxHairpin': 1,
            'loopMin': 6,
            'loopMax': 7,
            'homVisible': False,
            'hairpinVisible': False,
            'gcVisible': False,
            }

    if request.method == 'POST':

        if "analyseSeqSubmission" in request.form:
            # Form being submitted; grab data from form.
            gcContentMinPercentage = 25
            gcContentMaxPercentage = 65
            max_homopolymer = 1
            maxHairpin = 1
            loopMin = 6
            loopMax = 7
            constraints = set()
            if request.form['homVisible'] == 'True':
                constraints.add('hom')
                max_homopolymer = int(request.form['maxHomopolymer'])
            if request.form['hairpinVisible'] == 'True':
                constraints.add('hairpin')
                maxHairpin = int(request.form['maxHairpin'])
                loopMin = int(request.form['loopMin'])
                loopMax = int(request.form['loopMax'])
            if request.form['gcVisible'] == 'True':
                constraints.add('gcContent')
                gcContentMinPercentage = int(request.form['gcContentMinPercentage'])
                gcContentMaxPercentage = int(request.form['gcContentMaxPercentage'])

            form = {'gcContentMinPercentage': gcContentMinPercentage, 
                    'gcContentMaxPercentage': gcContentMaxPercentage, 
                    'maxHomopolymer': max_homopolymer,
                    'maxHairpin': maxHairpin,
                    'loopMin': loopMin,
                    'loopMax': loopMax,
                    'homVisible': request.form['homVisible'],
                    'hairpinVisible': request.form['hairpinVisible'],
                    'gcVisible': request.form['gcVisible'], 
                    }

            keys = request.form["keys"]
            payloads = request.form["payloads"]
            isValid = False

            message, keys, payloads, isValid = analyse_input.validate_keys_and_payloads(keys, \
                                                                    payloads, form, constraints)

            # Render the page
            return render_template('index.html', payloads=payloads, message=message, keys=keys, \
                                    form=form, isValid=isValid)
        
        return render_template('index.html', payloads="", keys="", message="", form=form, \
                                isValid=False)

    return render_template('index.html', payloads="", keys="", message="", form=form, \
                            isValid=False)


if __name__ == "__main__":
   app.run(debug=True)


