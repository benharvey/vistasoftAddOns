# VistasoftAddOns
Additions to Vistasoft MRI analysis repository to run quantity-selective neural response models

We have been using Vistasoft (https://github.com/vistalab/vistasoft) for fMRI analysis for years, and have extended it in several ways to allow the analyses underlying the studies listed below. 
Our main uses are to model 1-dimensional response functions (for example numerosity and size tuning), monotonic response models (for example monotonic responses to numerosity, Fourier power or event timing) or complex response models to events with variable timing.
But we have unfortunately found that our additions can affect the basic functionality of Vistasoft for visual field mapping (the most common use case), so it is better not to commit our changes to the main repository.
Therefore, this repository contains all the Vistasoft functions we have changed, and others we have added. If you want to run our non-visual-field-mapping response models, Vistasoft should be added to your Matlab path first, followed by this repository (VistasoftAddOns). This will then use the functions here in place of Vistasoft's regular functions.
Please note that computational modelling of fMRI data is very complex (even for visual field mapping, and particularly for timing-selective response models) and fMRI experiments should be designed from the start (before data collection) to suit this type of analysis.
Therefore, while we want to make this code public, it is a good idea to contact the authors (particularly Ben M Harvey) for advice about your experiment and analyses.

Some papers that have relied on the code in this repository:
Harvey BM, Klein BP, Petridou N, Dumoulin SO (2013) Topographic representation of numerosity in human parietal cortex. Science 341: 1123-1126.
Harvey BM, Fracasso A, Petridou N, Dumoulin SO (2015) Topographic representations of object size and relationships with numerosity reveal generalized quantity processing in human parietal cortex. Proceedings of the National Academy of Sciences 112(44): 13525-13530.
Harvey BM, Dumoulin SO (2016) Visual motion transforms visual space representations similarly throughout the human visual hierarchy. Neuroimage 127: 173-185.!
Harvey BM, Dumoulin SO (2017) Can responses to basic non-numerical visual features explain neural numerosity responses? Neuroimage 149: 200-209.
Harvey BM, Dumoulin SO (2017) A network of topographic numerosity maps in human association cortex. Nature Human Behaviour 1: 36.
Harvey BM, Dumoulin SO, Fracasso A, Paul JM (2020) A network of topographic maps in human association cortex hierarchically transforms visual timing-selective responses. Current Biology. 30: 1424-1434.
Hofstetter S, Cai Y, Harvey BM, Dumoulin SO (2021) Topographic maps representing haptic numerosity reveals distinct sensory representations in supramodal networks. Nature Communications. 12: 221.
Cai Y, Hofstetter S, van Dijk J, Zuiderbaan W, van der Zwaag W, Harvey BM, Dumoulin SO (2021) Numerosity maps cover subitizing and estimation ranges. Nature Communications, 12: 3374.
Tsouli A, Cai Y, van Ackooij M, Hofstetter S, Harvey BM, Te Pas S, van der Smagt, Dumoulin SO (2021) Adaptation to visual numerosity changes neural numerosity selectivity. NeuroImage, 229: 117794
Hendrikx E, Paul JM, van Ackooij M, van der Stoep N, Harvey BM (2021) Visual timing-tuned responses in human association cortices and response dynamics in early visual cortex. Nature Communications (In revision). Also uses: https://github.com/evi-hendrikx/MonoTunedTiming
Paul JM, van Ackooij M, ten Cate TC, Harvey BM (2021) Numerosity tuning in human association cortices and local image contrast representations in early visual cortex. Nature Communications (In press).
