/*
 * \file Chromosome.h
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

#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include <vector>
#include "gsl/gsl_rng.h"
using namespace std;

class Chromosome
{
  int nbSites;
  float probTEPerSite;
  int verbose;
  gsl_rng * r;

  vector<int> vSeq;

 public:
  Chromosome( void );
  Chromosome( int, float, int, gsl_rng* );
  bool operator==( const Chromosome & );
  Chromosome& operator=( const Chromosome& );
  int& operator[]( int );
  void reset( void );

  void setNbSites( int );
  void setProbTEsPerSite( float );
  void setVerbose( int );
  void setRng( gsl_rng * );
  void setSequence( vector<int> );

  int getNbSites( void );
  float getProbTEsPerSite( void );
  int getVerbose( void );
  gsl_rng* getRng( void );

  void initialize( void );
  int getNbTEs( void );
  void loss( void );
  void transposition( void );
  void printSequence( void );
  bool isTranspElemAtSite( int );
};

#endif
