# WPCA
# Robust weighted principal components analysis demodulation algorithm for phase-shifting interferometry

In this work, we improve our previous principal component analysis demodulation approach [1]. We present an asynchronous phase-shifting demodulation approach based on the principal component analysis demodulation method that is robust to typical problems as turbulence, vibrations, and temporal instabilities of the optical setup,  which represent the current main limitations of phase-shifting interfometry. The method brings together a two-step and a phase-shifting asynchronous demodulation method to share their benefits while reducing their intrinsic limitations. Thus, the proposed approach is based on a two-fold process. First, the modulating phase is estimated from a two-step demodulation approach. Second, this information is used to compute weights to each phase-shifted pattern of the interferogram sequence, which are used in a novel weighted principal component demodulation approach. The proposed technique has been tested with simulated and real interferograms affected by turbulence and vibrations providing very satisfactory results in challenging cases.

Here, we provide the code that reproduces all the results and figures shown in the paper [2]. This software has been tested under Matlab R2018a, R2019b, R2020b. Additionaly, it might be necessary to have instaled the Matlab Image Processing toolbox. To run the different examples, run the script WPCA_demo.m changing the variable 'medida' to values between 1 and 7. The cases medida=1,2,4 reproduce simulated results shown in [2], and medida=5,7 experimental cases.

The experimental images used in the script WPCA_demo.m and shown in the paper [2] can be downloaded from [3], [4] and [5].

[1] J Vargas, JA Quiroga, T Belenguer, "Phase-shifting interferometry based on principal component analysis", Optics Letters 36(8) 1326-1328 (2011)

[2] "Robust weighted principal components analysis demodulation algorithm for phase-shifting interferometry", (2020)

[3] https://ucomplutense-my.sharepoint.com/:u:/g/personal/jvargas_ucm_es/ERT0GRD7G5RClIjWJzk2OwAB3RBjW716hsgPirgGNy7Znw?e=kkqOMV

[4] https://ucomplutense-my.sharepoint.com/:u:/g/personal/jvargas_ucm_es/EdZrb39FG6hGmVr-HOtE8SwBuYI0CLUsXCM8P_xArV84Ow?e=gAiioV

[5] https://ucomplutense-my.sharepoint.com/:u:/g/personal/jvargas_ucm_es/EVOSW0Hsbc5FvEZzB5SIyHEBfF7UmLXAfkEi1voNq7ElGw?e=aUJqca 


