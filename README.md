# IVTmechanisticmodel

A mechanistic model for predicting reaction rates, RNA yields, and the presence of magnesium pyrophosphate crystallization in the in vitro transcription (IVT) reaction. While this is written in the Julia programming language, it is designed to be generally usable without a specific background in Julia or programming in general. The formulation and parameter estimation of the model is outlined in the publication:
Nathan Merica Stover, Krystian Ganko, Richard Braatz. Mechanistic Modeling of In Vitro Transcription. Authorea. August 25, 2023.

Specifically, this model simulates yields of the IVT reaction as a function of the following inputs:
1. T7 RNAP concentration
2. DNA template concentration
3. Concentration of each NTP (Including cap analogues)
4. Magnesium salt concentration
5. Inorganic pyrophosphatase concentration
6. Reaction time

In the notebook akamafittingdemo.ipynb, we demonstrate the methods used to generate the results and figures shown in our paper. However, this package includes a suite of tools that can be used by practioners to extend this model to their own work. This package can be used to:

1. Predict the results of new experiments
2. Visually compare model predictions and experimental data. This can be useful for gaining a mechanistic understanding  of why the IVT reaction behaves the way that it does. 
3. Fine-tune the model on new data. Specifically, this means performing the parameter estimation process performed in our origianl work on a new dataset comprising the original data and a new dataset. 

## Setting up package. 

This package is written in the Julia language, which you will need to download at https://julialang.org/downloads/ 
To use our demo files, you will also need to install jupyter notebooks: https://jupyter.org/
Next, download (or clone/fork if you are a pro) this repository to a desired location to your local machine. At this point you should be able to run the two ipynb files in the jupyter notebook interface. The first time you run should take a long time as it involves installing the required packages. 

## Analyzing new data with model

With this package, you can input your own experimental data for analysis. See data/pfizerNTP.csv as an example of the data input format. Each row represents a set of concentration inputs (i.e., an instance of running the IVT reaction.) The columns represent:

1. Species: An integer representing the chemical species that is being measured. Almost all data will be of RNA concentration measurements, which correponds to "2". However, the full list of possible inputs is
- DNA: 1
- RNA: 2
- Pyrophosphate: 3
- ATP: 4
- UTP: 5
- CTP: 6
- GTP: 7
- Mg: 8
- Phosphate: 10

2. Times: The time (in hours) of data collection. If only one time was collected, it can be written as a plain number. However, if multiple timepoints were collected for the same reaction set, of initial conditions (e.g. kinetics experiments), the times are written seperated by commas in parenthesis. For example, (0.5,1.0,2.0).

3. T7RNAP: The initial concentration of T7 RNA polymerase in nM. Note that we use the conversion factor 1 U/uL = 16 nM in our work.
4. DNA: Concentration of DNA template in nM
5. ATP/UTP/CTP/GTP/Cap: Initial concentration of respective nucleotide (or cap analogue) in mM
6. Mg: Initial concentration of Mg salt in mM
7. PPiase: Initial concentration of inorganic pyrophosphatase in (U/uL). Note that different suppliers define enzyme activity units differently. If your inorganic pyrophosphatase comes from Hongene, New England Bio, or any other company that defines U as "the amount of enzyme that will generate 1 µmol of phosphate per minute..", you need to divide your inorganic pyrophosphatase concentration by two to be congruent with the definition used in our model. 
8. Buffer: Concentration of buffer used (This is almost always 40 mM). Our model assumes that Tris-HCl buffer at pH ~8.0 is used. 
9. DTT/Spermadine: These components are not actually used in the model, but they appear in this sheet for ease of use in other analysis since they are often included in IVT schemes. You do not need to input real values into these cells. 
10. NA/NU/NC/NG: Number of each nucleotide in the RNA sequence of interest. 
11. output yield: Same format as times, but with the output measurement of interest. If the species measured is RNA, input in units of micromoles per liter. For all other species, input in units of millimoles per liter. 
12. outputstd: 1-sigma error bars of measurements. If you input a single number, it will be applied to all data in that row. If you have some sort of idea that some measurments are more accurate than others (this is rare), you can input using same format as the output yield column, where each item in the outputstd list is the error in the corresonding entry in the output yield list. DO NOT SET THIS VALUE TO BE ZERO, as it will introduce divide by zero errors.
"

## Fitting Model to New Data

You can use this model to analyze new data and test hypotheses (for example, how would yields change if I removed the pyrophosphatase enzyme?). Because our model is based in a first priciples understanding of the IVT reaciton, we expect that it is accurate about major trends (For example, the sensitivity of the reaction with respect to each of the components). However, because model prediction is unreliable when extrapolating to unexplored input regimes, the exact results of the model will commonly be incorrect when applied to new data. It is not fun to do hypothesis testing when the model does not describe the data in the first place! To help with this, we have included a framework to "fine-tune" the model on your custom data. Specifically, this adds your custom data to the origial calibration set and re-performs our parameter estimation process.

The notebook customfittingdemo.ipynb includes an example of this fine-tuning. Recent work from Pfizer (Aritra Sarkar, Guogang Dong, Jennifer Quaglia-Motta, Kelly Sackett, Flow-NMR as a process-monitoring tool for mRNA IVT reaction, Journal of Pharmaceutical Sciences,2023) has used novel NMR techniques to measure concentrations of each NTP and the phosphate ion with high sampling frequencies and low error. While the mechanistic model presented in our work was not calibrated using this kind of data, it can be incorporated in our mechanistic framework. 

We demonstrate how to fit our model to new data in the notebook customfittingdemo.jl. We first attempt to predict these new data using the parameters estimated in our original work. 

![pfizerdata_unfitted](https://github.com/nathanmstover/IVTmechanisticmodel/assets/97487659/0eb4472e-d878-4481-bf07-d96965b7998c)

Note that these predictions are incorrect with high uncertainty. Specifically, they are underpredicting the rate of reaction. The key issue here is that we are extrapolating in input space from the data used for calibration. Consider the Akama data of initial reaction rates that was used to develop these parameters. 

![p1](https://github.com/nathanmstover/IVTmechanisticmodel/assets/97487659/24d200d1-182c-4e4a-8a98-e35fce641ec4)


The reaction we are attempting to validate our model on used a total NTP concentration of about 25 mM and a Mg salt concentration of 16.5 mM. This is pretty far outside of the data range that we are calibrating our model on. Visually, it is clear that it is difficult to infer values at 25 mM NTP from the data we have. 

Now we consider adding the Pfizer data to our calibration data and rerunning our parameter estimation process. We get the following fits to the Pfizer data:

![pfizerdata_fitted](https://github.com/nathanmstover/IVTmechanisticmodel/assets/97487659/f899cc40-1d15-4cf4-a902-ed7ffa697a7e)

Pretty good! We can see how these new parameters affect other parts of our original calibration data. Our predictions for reaction rates are slightly revised.

![p2](https://github.com/nathanmstover/IVTmechanisticmodel/assets/97487659/6933f58e-c50e-495c-8a9b-d91708d175e8)

In general, we are fitting data collected at high (<10 mM) Mg concentrations better than before, and data at low Mg concentrations worse. The mechanism developed in our work appears to be unable describe both of these regimes accurately, which is the subject of future work. 

## Hypothesis Testing 
Now that we have a model that reasonably fits the entire data set we have thrown at it, we can consider using it to answer interesting questions. For example, we can test the effect of removing the in pyrophosphatase enzyme on the pfizer reaction conditions. We simply run the model on a new data file that does not include any real data. See data/pfizertest.csv as an example for how to we do this (we set the time column to the time window we are interested in studying this reaction for, and the output yield column to zero). Our comparison shows that without the pyrophosphatase enzyme, early reaction halting decreases RNA yield, which is consistent with the main thrust of our publication.

![ppiasehypothesistesting](https://github.com/nathanmstover/IVTmechanisticmodel/assets/97487659/dd3896b7-d339-48f2-9167-55db0938e694)

## Caveats and Best Practices with Model

When attempting to apply a mechanistic model like this to new data, there are a number of factors that can lead to incorrect results:
1. There is still a limited understanding of how sequence identity affects the transcription process in the context of biomanufacturing. Different mRNA sequences could transcribe at different rates, for example.
2. This model was developed using a dataset containing relatively low concentrations of NTPs and Mg. Extrapolating to new input spaces should not be expected to work perfectly.
3. There are a number of unknown unknowns when it comes to comparing the work of different researchers (analytical techniques, mixing, etc..)

In the face of all of this uncertainity, our recomendation of how best to use this tool is as follows:
1. First, attempt to predict your data using the original parameters estimated in our work. If the model does a good job of describing the data, you can consider using the model as a tool for hypothesis testing (in the local range of data inputs that you have - again, you should not expect the model to be accurate in extrapolation). 
2. If the model poorly predicts results or predicts with extremely high uncertainty, follow the fine tuning procedure.
3. If fine tuning results in good fits for the entirity of the dataset (your new data + the original data) you can consider using it as a platform for hypothesis testing (again, while avoiding extrapolation)



