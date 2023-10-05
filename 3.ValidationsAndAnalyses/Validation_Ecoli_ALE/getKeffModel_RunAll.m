%% Cleaning the workspace and the command window
clear;clc

%% Loading the enzyme-constrained model and general data
load('../../1.Minimal_adjustment/Code_Ecoli/eciML1515_batch.mat');

%% Generating keff adjusted model
current = pwd;
cd ../../1.Minimal_adjustment/PRESTO/Program

ecModel_batch.name = 'Escherichia coli';
enzMetPfx = 'prot_';

ecModel_batch_keff = create_pFBAmod(ecModel_batch, enzMetPfx);

cd (current)