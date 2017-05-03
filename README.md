# Redundant-Mutually-Connected-Giant-Component

 This code can be redistributed and/or modified
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at
  your option) any later version.
   
  This program is distributed ny the authors in the hope that it will be 
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 
   
  If you use this code please cite the following paper:
 
 [1] F. Radicchi and G. Bianconi, "Redundant interdependencies boost the robustness of multiplex networks" Phys. Rev. X 7, 011013 (2017)
 
 
  (c) G. Bianconi (email: ginestra.bianconi@gmail.com ) 
 
 ****************************************************************************************************************************************
 This program reads data (edgelist) for  a multiplex network with M=3 layers 
 in the format :
 layer node1 node2 weight (weight to be disregarded)
 
 It consider a single random configuration of initial damage  
 and evaluates the Redudant Mutually Connected Giant Component as a function of the fraction of  nodes damaged
 The Message Passing algorithm for each single realization of the random damage
 is calculated on the same network.
 
 INPUT
 N total number of nodes filename of the edgelist

 
 OUTPUT
 
 The output file "Redundant_MP.txt"  contains two columns: 
 p=1-f the fraction of non damaged nodes, and the size of the RMCGC predicted by the Message Passing Algorithm
 
 The output file "Redundant_sim.txt"  contains two columns: 
 p=1-f the fraction of non damaged nodes  and the size of the RMCGC obtained by direct numerical calculations
