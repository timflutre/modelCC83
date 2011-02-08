/*
 * \file Population.cpp
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
#include <iomanip>  // for setprecision
#include <fstream>
#include <numeric>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <typeinfo>  // for typeid
using namespace std;

#include "Population.h"
#include "Individual.h"

Population::Population( void )
{
  reset();
}

void Population::reset( void )
{
  setNbDiploids( 0 );
  setNbChrPerIndividual( 0 );
  setNbSitesPerChromosome( 0 );
  setExpNbTEsPerIndividual( 0 );
  setTotalMapDist( 0 );
  setZygoteSelection( false );
  setSelMultiplicator( 0.0 );
  setSelExponent( 0.0 );
  setVerbose( 0 );
  vInd.clear();
}

void Population::setNbDiploids( int nd )
{
  nbDiploids = nd;
}

void Population::setNbChrPerIndividual( int cpi )
{
  nbChrPerInd = cpi;
}

void Population::setNbSitesPerChromosome( int spc )
{
  nbSitesPerChr = spc;
}

void Population::setExpNbTEsPerIndividual( int nti )
{
  expNbTEsPerInd = nti;
}

void Population::setTotalMapDist( int tmd )
{
  totalMapDist = tmd;
}

void Population::setZygoteSelection( bool zs )
{
  zygoteSelection = zs;
}

void Population::setSelMultiplicator( float sm )
{
  selMult = sm;
}

void Population::setSelExponent( float se )
{
  selExp = se;
}

void Population::setVerbose( int v )
{
  verbose = v;
}

void Population::setRng( gsl_rng * rng )
{
  r = rng;
}

int Population::getNbDiploids( void )
{
  return( nbDiploids );
}

int Population::getNbChrPerIndividual( void )
{
  return( nbChrPerInd );
}

int Population::getNbSitesPerChromosome( void )
{
  return( nbSitesPerChr );
}

int Population::getExpNbTEsPerIndividual( void )
{
  return( expNbTEsPerInd );
}

int Population::getTotalMapDist( void )
{
  return( totalMapDist );
}

bool Population::getZygoteSelection( void )
{
  return( zygoteSelection );
}

float Population::getSelMultiplicator( void )
{
  return( selMult );
}

float Population::getSelExponent( void )
{
  return( selExp );
}

int Population::getVerbose( void )
{
  return( verbose );
}

gsl_rng* Population::getRng( void )
{
  return( r );
}

void Population::initialize( void )
{
  if( getVerbose() > 0 )
    cout << "initialization" << endl;
  for( int i=0; i<nbDiploids; ++i ){
    if( getVerbose() > 1 )
      cout << "initialize individual " << i+1 << endl;
    Individual ind;
    ind.setNbChromosomes( nbChrPerInd );
    ind.setNbSitesPerChromosome( nbSitesPerChr );
    ind.setExpNbTEsPerIndividual( expNbTEsPerInd );
    ind.setZygoteSelection( zygoteSelection );
    ind.setSelMultiplicator( selMult );
    ind.setSelExponent( selExp );
    ind.setVerbose( verbose-1 );
    ind.setRng( r );
    ind.initialize();
    vInd.push_back( ind );
  }
}

vector<double> Population::getNbTEsPerInd( void )
{
  vector<double> vNbTEsPerInd;
  for( int i=0; i<nbDiploids; ++i )
    vNbTEsPerInd.push_back( vInd[i].getNbTEs() );
  return( vNbTEsPerInd );
}

void Population::getNbTEsPerInd( gsl_vector *vNbTEsPerInd )
{
  for( int i=0; i<nbDiploids; ++i )
    gsl_vector_set( vNbTEsPerInd, i, vInd[i].getNbTEs() );
}

int Population::getSumNbTEs( void )
{
  vector<double> vNbTEsPerInd = getNbTEsPerInd();
  return( getSumNbTEs( vNbTEsPerInd ) );
}

int Population::getSumNbTEs( vector<double> vNbTEsPerInd )
{
  return( accumulate( vNbTEsPerInd.begin(),
                      vNbTEsPerInd.end(),
                      0 ) );
}

int Population::getSumNbTEs( gsl_vector_view gvNbTEsPerInd )
{
  return( gsl_blas_dasum( &gvNbTEsPerInd.vector ) );
}

float Population::getMeanNbTEs( gsl_vector_view gvNbTEsPerInd )
{
  return( gsl_stats_mean( gvNbTEsPerInd.vector.data, 1, nbDiploids ) );
}

float Population::getVarNbTEs( gsl_vector_view gvNbTEsPerInd )
{
  return( gsl_stats_variance( gvNbTEsPerInd.vector.data, 1, nbDiploids ) );
}

float Population::getSdNbTEs( gsl_vector_view gvNbTEsPerInd )
{
  if( getMeanNbTEs( gvNbTEsPerInd ) == 0 )
    return( 0 );
  return( gsl_stats_sd( gvNbTEsPerInd.vector.data, 1, nbDiploids ) );
}

int Population::getMinNbTEs( gsl_vector_view gvNbTEsPerInd )
{
  return( gsl_stats_min( gvNbTEsPerInd.vector.data, 1, nbDiploids ) );
}

float Population::getQuantileNbTEs( gsl_vector_view gvNbTEsPerInd, float q )
{
  gsl_sort_vector( &gvNbTEsPerInd.vector );
  return( gsl_stats_quantile_from_sorted_data( gvNbTEsPerInd.vector.data,
                                               1, nbDiploids, q ) );
}

int Population::getMaxNbTEs( gsl_vector_view gvNbTEsPerInd )
{
  return( gsl_stats_max( gvNbTEsPerInd.vector.data, 1, nbDiploids ) );
}

void Population::printDistribTEsPerInd( void )
{
  vector<double> vNbTEsPerInd = getNbTEsPerInd();
  gsl_vector_view gvNbTEsPerInd = gsl_vector_view_array( &vNbTEsPerInd[0],
                                                         vNbTEsPerInd.size() );
  cout << "TEs=" << getSumNbTEs( gvNbTEsPerInd );
  cout << " mean=" << setprecision(3) << getMeanNbTEs( gvNbTEsPerInd );
  cout << " sd="<< setprecision(3) << getSdNbTEs( gvNbTEsPerInd );
  cout << " min=" << getMinNbTEs( gvNbTEsPerInd );
  cout << " q25=" << setprecision(3) << getQuantileNbTEs( gvNbTEsPerInd, 0.25 );
  cout << " med=" << setprecision(3) << getQuantileNbTEs( gvNbTEsPerInd, 0.50 );
  cout << " q75=" << setprecision(3) << getQuantileNbTEs( gvNbTEsPerInd, 0.75 );
  cout << " max=" << getMaxNbTEs( gvNbTEsPerInd );
  cout << endl;
}

void Population::sampleCouple( Individual &parent1, Individual &parent2 )
{
  int idPar1 = gsl_rng_uniform_int( r, nbDiploids );
  int idPar2 = gsl_rng_uniform_int( r, nbDiploids );
  while( idPar2 == idPar1 )
    idPar2 = gsl_rng_uniform_int( r, nbDiploids );
  parent1 = vInd[ idPar1 ];
  parent2 = vInd[ idPar2 ];
}

void Population::addIndividual( void )
{

}

void Population::setIndividuals( vector<Individual> vNewInd )
{
  if( vNewInd.size() != (unsigned) getNbDiploids()
      || vNewInd[0].getNbChromosomes() != getNbChrPerIndividual()
      || vNewInd[0].getNbSitesPerChromosome() != getNbSitesPerChromosome() ){
    cerr << "ERROR: new population has different features" << endl;
    exit( EXIT_FAILURE );
  }
  vInd = vNewInd;
}

void Population::makeNewGeneration( int v )
{
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  vector<Individual> vNewInd;
  int i = 0;
  while( i < getNbDiploids() ){
    if( getVerbose() > 1 )
      cout << "make individual " << i+1 << endl << flush;
    Individual parent1, parent2;
    sampleCouple( parent1, parent2 );
    vector<Chromosome> gamete1, gamete2;
    parent1.getGamete( totalMapDist, gamete1 );
    parent2.getGamete( totalMapDist, gamete2 );
    Individual ind = Individual();
    ind.fecundation( gamete1, gamete2, zygoteSelection,
                     selMult, selExp, verbose-1 );
    if( ind.isViable() ){
      vNewInd.push_back( ind );
      ++i;
    }
  }
  setIndividuals( vNewInd );
}

void Population::loss( float probLoss )
{
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  int nbLosses = 0;
  for( int i=0; i<nbDiploids; ++i )
    nbLosses += vInd[i].loss( probLoss );
  if( getVerbose() > 0 )
    cout << "nb of losses: " << nbLosses << endl;
}

void Population::transposition( float probTransp0, float k )
{
  if( getVerbose() > 0 )
    cout << typeid(this).name() << "::" <<  __FUNCTION__ << endl << flush;
  int nbTransp = 0;
  for( int i=0; i<nbDiploids; ++i )
    nbTransp += vInd[i].transposition( probTransp0, k );
  if( getVerbose() > 0 )
    cout << "nb of transpositions: " << nbTransp << endl;
}

void Population::saveData( int simu, int gen, string outFile )
{
  string sep = "\t";
  ofstream outStream;
  outStream.open( outFile.c_str(),
                  fstream::in | fstream::out | fstream::app );
  outStream << simu << sep << gen << sep;

  vector<double> vNbTEsPerInd = getNbTEsPerInd();
  gsl_vector_view gvNbTEsPerInd = gsl_vector_view_array( &vNbTEsPerInd[0],
                                                         vNbTEsPerInd.size() );
  outStream << getSumNbTEs( gvNbTEsPerInd ) << sep;
  outStream << setprecision(3) << getMeanNbTEs( gvNbTEsPerInd ) << sep;
  outStream << setprecision(3) << getVarNbTEs( gvNbTEsPerInd ) << sep;
  outStream << setprecision(3) << getSdNbTEs( gvNbTEsPerInd ) << sep;
  outStream << getMinNbTEs( gvNbTEsPerInd ) << sep;
  outStream << setprecision(3) << getQuantileNbTEs( gvNbTEsPerInd, 0.25 ) << sep;
  outStream << setprecision(3) << getQuantileNbTEs( gvNbTEsPerInd, 0.50 ) << sep;
  outStream << setprecision(3) << getQuantileNbTEs( gvNbTEsPerInd, 0.75 ) << sep;
  outStream << getMaxNbTEs( gvNbTEsPerInd ) << sep;
  outStream << setprecision(3) << getPropEmptyLoci() << sep;
  outStream << endl;
  outStream.close();
}

void Population::getOccPerLocus( vector< vector<int> > & vOcc )
{
  for( int ind=0; ind<nbDiploids; ++ind )
    vInd[ ind ].getOccPerLocus( vOcc[ ind ] );
}

void Population::getFreqBetweenLoci( void )
{
  
}

float Population::getPropEmptyLoci( void )
{
  int nbLociPerInd = ( nbChrPerInd * nbSitesPerChr ) / 2;
  vector< vector<int> > vOcc ( nbDiploids, vector<int>( nbLociPerInd, 0 ) );
  getOccPerLocus( vOcc );
  int nbEmptyLoci = 0;
  for( int loc=0; loc<nbLociPerInd; ++loc )
    for( int ind=0; ind<nbDiploids; ++ind )
      if( vOcc[ ind ][ loc ] == 0 )
        ++ nbEmptyLoci;
  return( (float) nbEmptyLoci / ( nbLociPerInd * nbDiploids ) );
}
