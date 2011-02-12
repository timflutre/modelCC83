#include <iostream>
#include "gsl/gsl_rng.h"
using namespace std;

#include "Population.h"
#include "Individual.h"
#include "Chromosome.h"

void usage( char *program_name, int status )
{
  cerr << "usage: " << program_name << " [options]\n";
  cerr << "options:" << endl;
  cerr << "     -h: this help" << endl;
  cerr << "     -v: verbose (default=0/1/2)" << endl;
  exit( status );
}

int test_Individual_recombine( gsl_rng * r, int verbose )
{
  if( verbose > 0 ){
    cout << __FUNCTION__<< ": ";
    if( verbose > 1 )
      cout << endl;
  }

  int nbSitesPerChr = 10;
  Individual ind;
  ind.setNbChromosomes( 2 );
  ind.setNbSitesPerChromosome( nbSitesPerChr );
  ind.setRng( r );
  ind.initialize();

  Chromosome chr1( nbSitesPerChr, 0.1, 0, r );
  Chromosome chr2( nbSitesPerChr, 0, 0, r );
  for( int i=0; i<nbSitesPerChr; ++i )
    chr2[i] = 1;
  if( verbose > 1 ){
    cout << "initChr1: ";
    chr1.printSequence();
    cout << "initChr2: ";
    chr2.printSequence();
  }

  Chromosome expChr1( nbSitesPerChr, 0.1, 0, r );
  Chromosome expChr2( nbSitesPerChr, 0.1, 0, r );
  for( int i=0; i<3; ++i )
    expChr2[i] = 1;
  for( int i=3; i<6; ++i )
    expChr1[i] = 1;
  for( int i=6; i<10; ++i )
    expChr2[i] = 1;
  if( verbose > 1 ){
    cout << "expChr1: ";
    expChr1.printSequence();
    cout << "expChr2: ";
    expChr2.printSequence();
  }

  ind.setVerbose( 0 );
  ind.recombine( 1, chr1, chr2 );
  if( verbose > 1 ){
    cout << "obsChr1: ";
    chr1.printSequence();
    cout << "obsChr2: ";
    chr2.printSequence();
  }

  if( expChr1 == chr1 && expChr2 == chr2 ){
    if( verbose > 0 )
      cout << "TRUE" << endl;
    return( 0 );
  }
  else{
    if( verbose > 0 )
      cout << "FALSE" << endl;
    return( 1 );
  }
}

int test_Population_getFreqTEsPerLocus( gsl_rng * r, int verbose )
{
  if( verbose > 0 ){
    cout << __FUNCTION__<< ": ";
    if( verbose > 1 )
      cout << endl;
  }

  Population pop;
  pop.setNbDiploids( 2 );
  pop.setNbChrPerIndividual( 4 );
  pop.setNbSitesPerChromosome( 4 );
  pop.setExpNbTEsPerIndividual( 4 );
  pop.setRng( r );
  pop.initialize();
  if( verbose > 1 )
    pop.printChrSequencesPerInd();

  vector<double> vExp;
  vExp.push_back( 0.25 );
  vExp.push_back( 0.25 );
  vExp.push_back( 0.25 );
  vExp.push_back( 0.25 );
  vExp.push_back( 0.25 );
  vExp.push_back( 0.5 );
  vExp.push_back( 0.25 );
  vExp.push_back( 0.5 );

  vector<double> vObs = pop.getFreqTEsPerLocus();
  if( verbose > 1 ){
    cout << "vObs:";
    for( vector<double>::iterator it=vObs.begin(); it!=vObs.end(); ++it )
      cout << " " << *it;
    cout << endl;
  }

  if( vExp == vObs ){
    if( verbose > 0 )
      cout << "TRUE" << endl;
    return( 0 );
  }
  else{
    if( verbose > 0 )
      cout << "FALSE" << endl;
    return( 1 );
  }
}

int test_Individual_getOccPerLocus( gsl_rng * r, int verbose )
{
  if( verbose > 0 ){
    cout << __FUNCTION__<< ": ";
    if( verbose > 1 )
      cout << endl;
  }

  Individual ind;
  ind.setNbChromosomes( 4 );
  ind.setNbSitesPerChromosome( 2 );
  ind.setExpNbTEsPerIndividual( 4 );
  ind.setRng( r );
  ind.initialize();
  if( verbose > 1 )
    ind.printChromosomes();

  int nbLoci = ind.getNbLoci();
  vector<int> vExp;
  vExp.push_back( 2 );
  vExp.push_back( 0 );
  vExp.push_back( 0 );
  vExp.push_back( 1 );

  vector<int> vObs( nbLoci, 0 );
  ind.getOccPerLocus( vObs );

  if( vExp == vObs ){
    if( verbose > 0 )
      cout << "TRUE" << endl;
    return( 0 );
  }
  else{
    if( verbose > 0 )
      cout << "FALSE" << endl;
    return( 1 );
  }
}

int test_Individual_getNbTEsForLocus( gsl_rng * r, int verbose )
{
  if( verbose > 0 ){
    cout << __FUNCTION__<< ": ";
    if( verbose > 1 )
      cout << endl;
  }

  Individual ind;
  ind.setNbChromosomes( 4 );
  ind.setNbSitesPerChromosome( 2 );
  ind.setExpNbTEsPerIndividual( 7 );
  ind.setRng( r );
  ind.initialize();
  if( verbose > 1 )
    ind.printChromosomes();

  int locus = 3;
  int exp = 1;

  int obs = ind.getNbTEsForLocus( locus );
  if( verbose > 1 )
    cout << "locus=" << locus << " nbTEsExp=" << exp << " nbTEsObs=" << obs << endl;

  if( exp == obs ){
    if( verbose > 0 )
      cout << "TRUE" << endl;
    return( 0 );
  }
  else{
    if( verbose > 0 )
      cout << "FALSE" << endl;
    return( 1 );
  }
}

int test_Individual_getNbSites( gsl_rng * r, int verbose )
{
  if( verbose > 0 ){
    cout << __FUNCTION__<< ": ";
    if( verbose > 1 )
      cout << endl;
  }

  Individual ind;
  ind.setNbChromosomes( 4 );
  ind.setNbSitesPerChromosome( 2 );
  ind.setExpNbTEsPerIndividual( 4 );
  ind.setRng( r );
  ind.initialize();
  if( verbose > 1 )
    ind.printChromosomes();

  int exp = 8;

  int obs = ind.getNbSites();
  if( verbose > 1 )
    cout << "nbSitesExp=" << exp << " nbSitesObs=" << obs << endl;

  if( exp == obs ){
    if( verbose > 0 )
      cout << "TRUE" << endl;
    return( 0 );
  }
  else{
    if( verbose > 0 )
      cout << "FALSE" << endl;
    return( 1 );
  }
}

int main( int argc, char* argv[] )
{
  int nbFalses = 0;
  int seed = 1859;
  int verbose = 0;
  int nbTests = 5;

  char c;
  extern char *optarg;
  while( (c = getopt(argc,argv,"hv:")) != -1 ){
    switch (c){
    case 'h':
      usage( argv[0], EXIT_SUCCESS );
      break;
    case 'v':
      verbose = atoi(optarg);
      break;
    case '?':
      usage( argv[0], EXIT_FAILURE );
    default:
      usage( argv[0], EXIT_FAILURE );
    }
  }

  gsl_rng * r;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc( T );
  gsl_rng_set( r, seed );

  nbFalses += test_Individual_recombine( r, verbose );
  nbFalses += test_Population_getFreqTEsPerLocus( r, verbose );
  nbFalses += test_Individual_getOccPerLocus( r, verbose );
  nbFalses += test_Individual_getNbTEsForLocus( r, verbose );
  nbFalses += test_Individual_getNbSites( r, verbose );

  cout << "errors: " << nbFalses
       << " / " << nbTests << endl;

  gsl_rng_free( r );
}
