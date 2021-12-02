============================= Assignment 5 - Digit Recognition using HMM =============================

Submission by : 214101011

To run : just load the project . build and run.

Folders in the project:
->digits : Holds the 10 digits to train , 20 utterance each
->tests : Holds the 10 digits of 10 utterances to test.
->models : holds the converged and averaged models of 10 digits.

Files in the project:
-HMM.cpp : building model and opimization.
-get_ceps.h : preprocess the samples and get universe of Cis
-lbg.h : from universe , generate codebook. 
-utility.h : some utilities.
-codebook.txt : codebook provided by TA. (Used this for this assignment , but built own codebook using LBG for final project. )
-Recording _Module.exe : for live testing.

Files generated during execution :
- output.txt : Shows all averaged models.
- universe.csv - Cis of all training digits
- Obs seq.csv - Final observation sequence from all training digits
- livetest.txt : Sample file generated while live testing.



