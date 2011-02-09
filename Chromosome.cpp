/*
 * \file Chromosome.cpp
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

#include <iostream>
#include <numeric>  // for accumulate
using namespace std;

#include "Chromosome.h"

Chromosome::Chromosome( void )
{
  reset();
}

Chromosome::Chromosome( int ns, float pts, int v, gsl_rng * rng )
{
  reset();
  nbSites = ns;
  probTEPerSite = pts;
  verbose = v;
  r = rng;
  vSeq.assign( nbSites, 0 );
}

bool Chromosome::operator==( const Chromosome &other )
{
  return( nbSites == other.nbSites
          && probTEPerSite == probTEPerSite
          && verbose == other.verbose
          && r == other.r
          && vSeq == other.vSeq );
}

Chromosome& Chromosome::operator=( Chromosome const& chr )
{
  nbSites = chr.nbSites;
  probTEPerSite = chr.probTEPerSite;
  verbose = chr.verbose;
  r = chr.r;
  vSeq = chr.vSeq;
  return( *this );
}

int& Chromosome::operator[]( int i )
{
  return( vSeq[i] );
}

void Chromosome::reset( void )
{
  setNbSites( 0 );
  setProbTEsPerSite( 0.0 );
  setVerbose( 0 );
  vSeq.clear();
}

void Chromosome::setNbSites( int ns )
{
  nbSites = ns;
}

void Chromosome::setProbTEsPerSite( float pts )
{
  probTEPerSite = pts;
}

void Chromosome::setVerbose( int v )
{
  verbose = v;
}

void Chromosome::setRng( gsl_rng * rng )
{
  r = rng;
}

void Chromosome::setSequence( vector<int> v )
{
  if( v.size() != (unsigned) nbSites ){
    cerr << "ERROR: try to initialize chromosome with bad sequence" << endl;
    exit( EXIT_FAILURE );
  }
  vSeq.clear();
  vSeq = v;
}

int Chromosome::getNbSites( void )
{
  return( nbSites );
}

float Chromosome::getProbTEsPerSite( void )
{
  return( probTEPerSite );
}

int Chromosome::getVerbose( void )
{
  return( verbose );
}

gsl_rng* Chromosome::getRng( void )
{
  return( r );
}

void Chromosome::initialize( void )
{
  vSeq.resize( nbSites );
  for( int i=0; i<nbSites; ++i ){
    if( getVerbose() > 0 )
      cout << "initialize site " << i+1 << endl;
    float probTE = gsl_rng_uniform( r );
    if( probTE < probTEPerSite )
      vSeq[i] = 1;
  }
}

int Chromosome::getNbTEs( void )
{
  return( accumulate( vSeq.begin(),
                      vSeq.end(),
                      0 ) );
}

void Chromosome::loss( void )
{
  int rankLostTE = gsl_rng_uniform_int( r, getNbTEs() ) + 1;
  int site = 0;
  int rankCurrentTE = 1;
  while( true ){
    if( vSeq[ site ] == 1 ){
      if( rankCurrentTE == rankLostTE ){
        vSeq[ site ] = 0;
        break;
      }
      else
        ++ rankCurrentTE;
    }
    ++ site;
  }
}

void Chromosome::transposition( void )
{
  int insSite = gsl_rng_uniform_int( r, nbSites );
  while( vSeq[ insSite ] == 1 )
    insSite = gsl_rng_uniform_int( r, nbSites );
  vSeq[ insSite ] = 1;
}

void Chromosome::printSequence( void )
{
  for( int i=0; i<nbSites; ++i )
    cout << vSeq[i];
  cout << endl;
}

bool Chromosome::isTranspElemAtSite( int site )
{
  return( vSeq[ site ] == 1 );
}
