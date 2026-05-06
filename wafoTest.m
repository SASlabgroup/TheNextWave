clear; close all; clc; 

addpath(genpath('C:\Users\Brittany Lydon\MREL Dropbox\Brittany Nicole Lydon\MATLAB\GitHub\wafo'))

load('.\ExampleData1\wavespec.mat')

S = struct();
S.Alt = wavespec.Etheta;
S.w = 2*pi/wavespec.f;
S.theta = wavespec.theta;

W = spec2field(S);