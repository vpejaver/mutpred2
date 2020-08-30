%
% function [prediction] = VLXT (sequence)
%
% Funtion predicts intrinsically disordered residues using VLXT classification model
%
% Input arguments:
%       sequence = string of amino acids
%                  sequence should be 30 symbols or more and should contain 20 amino acid codes only
%
% Output arguments:
%        prediction = vector of floats, approximates probabilities that each residue is disordered,
%                     residues with predicton above 0.5 are considered disordered
%
% Comment: in case of illegal characters, too short an input sequence etc. the program outputs error
%
% Reference for the VLXT predictor:
% P. Romero, Z. Obradovic, X. Li, E. C. Garner, C. J. Brown, A. K. Dunker (2001).
% Sequence complexity of disordered proteiin. Proteins: Structure, Function and Genetics, 42: 38-48.
%
% CPP function developed by: 
% Pedro Romero, 2000
% 
% Matlab function developed by: 
% Predrag Radivojac
% Indiana University
% Indianapolis, Indiana, May 2004
