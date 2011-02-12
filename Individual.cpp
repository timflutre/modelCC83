/*
 * \file Individual.cpp
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
#include <cmath>
#include <gsl/gsl_randist.h>
#include <typeinfo>
using namespace std;

#include "Individual.h"
#include "Chromosome.h"

Individual::Individual( void )
{
  reset();
}

Individual& Individual::operator=( Individual const& ind )
{
  nbChr = ind.nbChr;
  nbSitesPerChr = ind.nbSitesPerChr;
  expNbTEsPerInd = ind.expNbTEsPerInd;
  selMult = ind.selMult;
  selExp = ind.selExp;
  verbose = ind.verbose;
  r = ind.r;
  vChr = ind.vChr;
  return( *this );
}

void Individual::reset( void )
{
  setNbChromosomes( 0 );
  setNbSitesPerChromosome( 0 );
  setExpNbTEsPerIndividual( 0 );
  setSelMultiplicator( 0.0 );
  setSelExponent( 0.0 );
  setVerbose( 0 );
  vChr.clear();
}

void Individual::setNbChromosomes( int nc )
{
  nbChr = nc;
}

void Individual::setNbSitesPerChromosome( int spc )
{
  nbSitesPerChr = spc;
}

void Individual::setExpNbTEsPerIndividual( int nti )
{
  expNbTEsPerInd = nti;
}

void Individual::setZygoteSelection( bool zs )
{
  zygoteSelection = zs;
}

void Individual::setSelMultiplicator( float sm )
{
  selMult = sm;
}

void Individual::setSelExponent( float se )
{
  selExp = se;
}

void Individual::setVerbose( int v )
{
  verbose = v;
}

void Individual::setRng( gsl_rng * rng )
{
  r = rng;
}

void Individual::setChromosomes( vector<Chromosome> v )
{
  if( v.size() != (unsigned) nbChr ){
    cerr << "ERROR: different sizes in Individual::setChromosomes()" << endl;
    exit( EXIT_FAILURE );
  }
  vChr = v;
}

int Individual::getNbChromosomes( void )
{
  return( nbChr );
}

int Individual::getNbSitesPerChromosome( void )
{
  return( nbSitesPerChr );
}

int Individual::getExpNbTEsPerIndividual( void )
{
  return( expNbTEsPerInd );
}

float Individual::getSelMultiplicator( void )
{
  return( selMult );
}

float Individual::getSelExponent( void )
{
  return( selExp );
}

int Individual::getVerbose( void )
{
  return( verbose );
}

gsl_rng* Individual::getRng( void )
{
  return( r );
}

void Individual::initialize( void )
{
  float probTEPerSite = expNbTEsPerInd / float( nbChr * nbSitesPerChr );
  for( int i=0; i<nbChr; ++i ){
    if( getVerbose() > 0 )
      cout << "initialize chromosome " << i+1 << endl;
    Chromosome chr;
    chr.setNbSites( nbSitesPerChr );
    chr.setProbTEsPerSite( probTEPerSite );
    chr.setRng( r );
    chr.setVerbose( verbose-1 );
    chr.initialize();
    vChr.push_back( chr );
  }
}

int Individual::getNbTEs( void )
{
  int nbTEs=0;
  for( int i=0; i<nbChr; ++i )
    nbTEs += vChr[i].getNbTEs();
  return( nbTEs );
}

void Individual::getGamete( int totalMapDist, vector<Chromosome> &gamete )
{
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  recombine( totalMapDist, vChr[0], vChr[1] );
  int idChr1 = gsl_rng_uniform_int( r, 2 );  // chr from the 1st pair of homologues
  recombine( totalMapDist, vChr[2], vChr[3] );
  int idChr2 = gsl_rng_uniform_int( r, 2 ) + 2;  // chr from the 2nd pair of homologues
  gamete.clear();
  gamete.push_back( vChr[ idChr1 ] );
  gamete.push_back( vChr[ idChr2 ] );
}

void Individual::recombine( int totalMapDist, Chromosome &vChrA, Chromosome &vChrB )
{
  if( getVerbose() > 1 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  int nbCrossOvers = gsl_ran_poisson( r, totalMapDist );
  if( nbCrossOvers > 0 ){
    if( getVerbose() > 2 )
      cout << "nb of crossing-overs: " << nbCrossOvers << endl;
    for( int i=0; i<nbCrossOvers; ++i ){
      int coLocus = gsl_rng_uniform_int( r, nbSitesPerChr );
      if( getVerbose() > 3 ){
        cout << "crossing-over locus: " << coLocus+1 << endl;
        cout << "before crossing-over:" << endl;
        vChrA.printSequence();
        vChrB.printSequence();
      }
      Chromosome tmp = vChrA;
      for( int j=coLocus; j<nbSitesPerChr; ++j ){
        vChrA[j] = vChrB[j];
        vChrB[j] = tmp[j];
      }
      if( getVerbose() > 3 ){
        cout << "after crossing-over:" << endl;
        vChrA.printSequence();
        vChrB.printSequence();
      }
    }
  }
}

void Individual::fecundation( vector<Chromosome> gam1,
                              vector<Chromosome> gam2,
                              bool zs,
                              float sm,
                              float se,
                              int v )
{
  setVerbose( v );
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  setNbChromosomes( gam1.size() + gam2.size() );
  setNbSitesPerChromosome( gam1[0].getNbSites() );
  setZygoteSelection( zs );
  setSelMultiplicator( sm );
  setSelExponent( se );
  setRng( gam1[0].getRng() );
  vChr.push_back( gam1[0] );
  vChr.push_back( gam2[0] );
  vChr.push_back( gam1[1] );
  vChr.push_back( gam2[1] );
}

int Individual::loss( float probLoss )
{
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  int nbTEs = getNbTEs();
  if( nbTEs > 0 ){
    float meanNbLoss = probLoss * nbTEs;
    int nbLoss = gsl_ran_poisson( r, meanNbLoss );
    if( nbLoss > 0 ){
      if( getVerbose() > 1 )
        cout << "nb of losses: " << nbLoss << endl;
      for( int loss=0; loss<nbLoss; ++loss ){
        int chr = gsl_rng_uniform_int( r, nbChr );
        while( vChr[ chr ].getNbTEs() == 0 )
          chr = gsl_rng_uniform_int( r, nbChr );
        vChr[ chr ].loss();
      }
      if( getNbTEs() != nbTEs - nbLoss ){
        cerr << "ERROR: bad number of lost TEs (" << getNbTEs()
             << "!=" << nbTEs-nbLoss << ")" << endl;
        exit( EXIT_FAILURE );
      }
    }
    return( nbLoss );
  }
  return( 0 );
}

int Individual::transposition( float probTransp0, float k )
{
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  int nbTEs = getNbTEs();
  if( nbTEs > 0 ){
    float probTransp;
    if( k == 0 )
      probTransp = probTransp0;
    else
      probTransp = probTransp0 / float( 1 + k * nbTEs );
    float meanNbTransp = probTransp * nbTEs;
    int nbTransp = gsl_ran_poisson( r, meanNbTransp );
    if( nbTEs + nbTransp >= nbChr * nbSitesPerChr ){
      cerr << "WARNING: too many TEs and no more empty sites" << endl;
//       nbTransp = nbChr * nbSitesPerChr - nbTEs;
      exit( EXIT_FAILURE );
    }
    if( nbTransp > 0 ){
      if( getVerbose() > 1 )
        cout << "nb of transpositions: " << nbTransp << endl;
      for( int transp=0; transp<nbTransp; ++transp ){
        int chr = gsl_rng_uniform_int( r, nbChr );
        while( vChr[ chr ].getNbTEs() == nbSitesPerChr )
          chr = gsl_rng_uniform_int( r, nbChr );
        vChr[ chr ].transposition();
      }
      if( getNbTEs() != nbTEs + nbTransp ){
        cerr << "ERROR: bad number of transposed TEs (" << getNbTEs()
             << "!=" << nbTEs+nbTransp << ")" << endl;
        exit( EXIT_FAILURE );
      }
    }
    return( nbTransp );
  }
  return( 0 );
}

void Individual::getOccPerLocus( vector<int> & vOccInd )
{
  int locus = 0;
  for( int chr=0; chr<nbChr; chr+=2 )  // "+=2" -> diploids
    for( int site=0; site<nbSitesPerChr; ++site ){
      if( vChr[ chr ].isTranspElemAtSite( site ) )
        ++ vOccInd[ locus ];
      if( vChr[ chr+1 ].isTranspElemAtSite( site ) )
        ++ vOccInd[ locus ];
      ++ locus;
    }
}

float Individual::getFitness( void )
{
  return( 1 - selMult * pow( getNbTEs(), selExp ) );
}

bool Individual::isViable( void )
{
  if( not zygoteSelection )
    return( true );
  else{
    float probSel = gsl_rng_uniform( r );
    if( probSel <= getFitness() )
      return( true );
    else
      return( false );
  }
}

void Individual::printChromosomes( void )
{
  cout << "chromosomes (" << nbChr/2 << " pairs):" << endl;
  for( int i=0; i<nbChr; ++i )
    vChr[i].printSequence();
}

Chromosome& Individual::getChromosome( int idChr )
{
  return( vChr[idChr] );
}

int Individual::getNbTEsForLocus( int locus )
{
  int nbTEs = 0;
  int chrPair = floor( (float) locus / nbSitesPerChr );
  if( vChr[ 2*chrPair ].isTranspElemAtSite( locus-nbSitesPerChr*chrPair ) )
    ++ nbTEs;
  if( vChr[ 2*chrPair + 1 ].isTranspElemAtSite( locus-nbSitesPerChr*chrPair ) )
    ++ nbTEs;
  return( nbTEs );
}

int Individual::getNbLoci( void )
{
  if( vChr.size() == 0 )
    return( 0 );
  else if( vChr[0].getNbSites() == 0 )
    return( 0 );
  return( ( vChr.size() * vChr[0].getNbSites() ) / 2 );  // "/2" -> diploids
}

int Individual::getNbSites( void )
{
  if( vChr.size() == 0 )
    return( 0 );
  else if( vChr[0].getNbSites() == 0 )
    return( 0 );
  int nbSites = 0;
  for( int chr=0; chr<nbChr; ++chr )
    nbSites += vChr[ chr ].getNbSites();
  return( nbSites );
}
