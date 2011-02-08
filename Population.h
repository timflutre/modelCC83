/*
 * \file Population.h
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

#ifndef POPULATION_H
#define POPULATION_H

#include <vector>
#include <string>
#include <gsl/gsl_vector.h>
#include "gsl/gsl_rng.h"
using namespace std;

#include "Individual.h"

class Population
{
  int nbDiploids;
  int nbChrPerInd;
  int nbSitesPerChr;
  int expNbTEsPerInd;
  int totalMapDist;
  bool zygoteSelection;
  float selMult;
  float selExp;
  int verbose;
  gsl_rng * r;

  vector<Individual> vInd;

 public:
  Population( void );
  void reset( void );

  void setNbDiploids( int );
  void setNbChrPerIndividual( int );
  void setNbSitesPerChromosome( int );
  void setExpNbTEsPerIndividual( int );
  void setTotalMapDist( int );
  void setZygoteSelection( bool );
  void setSelMultiplicator( float );
  void setSelExponent( float );
  void setVerbose( int );
  void setRng( gsl_rng * );

  int getNbDiploids( void );
  int getNbChrPerIndividual( void );
  int getNbSitesPerChromosome( void );
  int getExpNbTEsPerIndividual( void );
  int getTotalMapDist( void );
  bool getZygoteSelection( void );
  float getSelMultiplicator( void );
  float getSelExponent( void );
  int getVerbose( void );
  gsl_rng* getRng( void );

  void initialize( void );
  vector<double> getNbTEsPerInd( void );
  void getNbTEsPerInd( gsl_vector * );
  int getSumNbTEs( void );
  int getSumNbTEs( vector<double> );
  int getSumNbTEs( gsl_vector_view );
  float getMeanNbTEs( gsl_vector_view );
  float getVarNbTEs( gsl_vector_view );
  float getSdNbTEs( gsl_vector_view );
  int getMinNbTEs( gsl_vector_view );
  float getQuantileNbTEs( gsl_vector_view, float );
  int getMaxNbTEs( gsl_vector_view );
  void printDistribTEsPerInd( void );
  void sampleCouple( Individual &, Individual & );
  void addIndividual( void );
  void setIndividuals( vector<Individual> );
  void makeNewGeneration( int );
  void loss( float );
  void transposition( float, float );
  void saveData( int, int, string );
  void getOccPerLocus( vector< vector<int> > & );
  void getFreqBetweenLoci( void );
  float getPropEmptyLoci( void );
};

#endif
