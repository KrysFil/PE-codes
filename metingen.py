from SpecClass import spectrum_analyzer
import transfers as tf

# Initialize the class, set up the models and change the parameters.
spek = spectrum_analyzer(2, 50000, 20, 5, 10, 10000, 2)
functionOfModel = None
model = tf.model(functionOfModel, R=None, C=None, L=None)
fitModel = tf.fitModel(functionOfModel)

# Measure and analyse with spectrum_analyzer methods
spek.meter("sessie6_1")                     # Measure the data
spek.magnitude_plotter(theory=False)        # Run the analysis, help me for difficult plots, set a freq from where begin
"Then"
# params,_ = spek.fitter(fitModel,[None,None])# Fit the data to the model
"Then"
# spek.model = model                          # Update the model manually after fit
# spek.magnitude_plotter(calcSlope=False)     # Plot a nice Bode plot now


