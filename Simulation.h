/*
 * \file Simulation.h
 */

// Purpose: simulate transposable elements dynamics in genomes with the 
// model of Charlesworth & Charlesworth (1983).
// Copyright (C) <2011>  <Timothée Flutre>
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include "gsl/gsl_rng.h"
using namespace std;

class Simulation
{
  int simuId;
  int nbGen;
  int nbDiploids;
  int nbChrPerInd;
  int nbSitesPerChr;
  int expNbTEsPerInd;
  int totalMapDist;
  float probLoss;
  float probTransp0;
  float k;
  bool zygoteSelection;
  float selMult;
  float selExp;
  string outFile;
  int verbose;
  gsl_rng * r;
  
 public:
  Simulation( void );

  void setSimulationIdentifier( int );
  void setNbGenerations( int );
  void setNbDiploids( int );
  void setNbChrPerIndividuals( int );
  void setNbSitesPerChromosome( int );
  void setExpNbTEsPerIndividual( int );
  void setTotalMapDist( int );
  void setProbLoss( float );
  void setProbTransp0( float );
  void setK( float );
  void setZygoteSelection( bool );
  void setSelMultiplicator( float );
  void setSelExponent( float );
  void setSeed( int );
  void setOutFile( string );
  void setVerbose( int );
  void setRng( gsl_rng * );

  int getSimulationIdentifier( void );
  int getNbGenerations( void );
  int getNbDiploids( void );
  int getNbChrPerIndividual( void );
  int getNbSitesPerChromosome( void );
  int getExpNbTEsPerIndividual( void );
  int getTotalMapDist( void );
  float getProbLoss( void );
  float getProbTransp0( void );
  float getK( void );
  bool getZygoteSelection( void );
  float getSelMultiplicator( void );
  float getSelExponent( void );
  string getOutFile( void );
  int getVerbose( void );

  void printSimGen( int );
  void run( void );
};

#endif
