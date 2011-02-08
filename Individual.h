/*
 * \file Individual.h
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

#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include <vector>
#include "gsl/gsl_rng.h"
using namespace std;

#include "Chromosome.h"

class Individual
{
  int nbChr;
  int nbSitesPerChr;
  int expNbTEsPerInd;
  bool zygoteSelection;
  float selMult;  // "s" parameter in Charlesworth & Charlesworth
  float selExp;   // "t" parameter in Charlesworth & Charlesworth
  int verbose;
  gsl_rng * r;

  vector<Chromosome> vChr;

 public:
  Individual( void );
  Individual& operator=( const Individual& );
  void reset( void );
  
  void setNbChromosomes( int );
  void setNbSitesPerChromosome( int );
  void setExpNbTEsPerIndividual( int );
  void setZygoteSelection( bool );
  void setSelMultiplicator( float );
  void setSelExponent( float );
  void setVerbose( int );
  void setRng( gsl_rng * );
  void setChromosomes( vector<Chromosome> );

  int getNbChromosomes( void );
  int getNbSitesPerChromosome( void );
  int getExpNbTEsPerIndividual( void );
  bool getZygoteSelection( void );
  float getSelMultiplicator( void );
  float getSelExponent( void );
  int getVerbose( void );
  gsl_rng* getRng( void );

  void initialize( void );
  int getNbTEs( void );
  void getGamete( int totalMapDist, vector<Chromosome> & );
  void recombine( int, Chromosome &, Chromosome & );
  void fecundation( vector<Chromosome>, vector<Chromosome>,
                    bool, float, float, int );
  int loss( float );
  int transposition( float, float );
  void getOccPerLocus( vector<int> & );
  float getFitness( void );
  bool isViable( void );
  void printChromosomes( void );
  Chromosome& getChromosome( int );
  int getNbTEsForLocus( int );
};

#endif
