/*
 * \file dynamics-TEs_cc-83.cpp
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
#include <iomanip>
#include <fstream>
#include <cstdlib>  // for EXIT_SUCCESS
#include <cstdio>  // for EOF
#include <sys/stat.h>  // for struct stat
#include <ctime>
#include "gsl/gsl_rng.h"
using namespace std;

#include "Simulation.h"

void usage( char *program_name, int status )
{
  cerr << "usage: " << program_name << " [options]\n";
  cerr << "options:" << endl;
  cerr << "     -h: this help" << endl;
  cerr << "     -s: number of simulations (default=1)" << endl;
  cerr << "     -n: number of diploids (default=10)" << endl;
  cerr << "     -g: number of generations per simulation (default=10)" << endl;
  cerr << "     -c: number of sites per chromosome (default=31)" << endl;
  cerr << "     -i: initial number of TEs per individual (default=10)" << endl;
  cerr << "     -t: transposition probability per TE per generation (default=0.01)" << endl;
  cerr << "     -k: parameter for transposition regulation (default=0.05)" << endl;
  cerr << "     -l: loss probability per TE per generation (default=0.005)" << endl;
  cerr << "     -d: total recombination map distance (default=90)" << endl;
  cerr << "         loose linkage: 90 units" << endl;
  cerr << "         tight linkage: 9 units" << endl;
  cerr << "     -S: apply zygote selection (eventually put k=0)" << endl;
  cerr << "     -m: selection multiplicator (only with -S, default=0.001)" << endl;
  cerr << "     -e: selection exponent (only with -S, default=1.5)" << endl;
  cerr << "     -r: seed of the pseudo-random generator (default=1859)" << endl;
  cerr << "     -o: name of the output file (default=data.csv)" << endl;
  cerr << "     -v: verbose (default=0/1/2)" << endl;
  exit( status );
}

void parse_args
( int argc, char **argv,
  int & nbSimu,
  int & nbDiploids,
  int & nbGen,
  int & nbSitesPerChr,
  int & initNbTEsPerInd,
  float & probTransp0,
  float & k,
  float & probLoss,
  int & totalMapDist,
  bool & zygoteSelection,
  float & selMult,
  float & selExp,
  int & seed,
  string & outFile,
  int & verbose
  )
{
  char c;
  extern char *optarg;
  while( (c = getopt(argc,argv,"hs:n:g:c:i:t:k:l:d:Sm:e:r:o:v:")) != -1 ){
    switch (c){
    case 'h':
      usage( argv[0], EXIT_SUCCESS );
      break;
    case 's':
      nbSimu = atoi(optarg);
      if( nbSimu <= 0 ){
        cerr << "ERROR: requires at least 1 simulation (-s)" << endl;
        usage( argv[0], EXIT_FAILURE );
      }
      break;
    case 'n':
      nbDiploids = atoi(optarg);
      if( nbDiploids <= 1 ){
        cerr << "ERROR: requires at least 2 individuals (-n)" << endl;
        usage( argv[0], EXIT_FAILURE );
      }
      break;
    case 'g':
      nbGen = atoi(optarg);
      break;
    case 'c':
      nbSitesPerChr = atoi(optarg);
      if( nbSitesPerChr <= 3 ){
        cerr << "ERROR: requires at least 3 sites per chromosome (-c)" << endl;
        usage( argv[0], EXIT_FAILURE );
      }
      break;
    case 'i':
      initNbTEsPerInd = atoi(optarg);
      if( initNbTEsPerInd <= 0){
        cerr << "ERROR: requires at least 1 TE (-i)" << endl;
        usage( argv[0], EXIT_FAILURE );
      }
      break;
    case 't':
      probTransp0 = atof(optarg);
      if( probTransp0 < 0 || probTransp0 > 1 ){
        cerr << "ERROR: probability should be between 0 and 1 (-t)" << endl;
        usage( argv[0], EXIT_FAILURE );
      }
      break;
    case 'k':
      k = atof(optarg);
      break;
    case 'l':
      probLoss = atof(optarg);
      if( probTransp0 < 0 || probTransp0 > 1 ){
        cerr << "ERROR: probability should be between 0 and 1 (-t)" << endl;
        usage( argv[0], EXIT_FAILURE );
      }
      break;
    case 'd':
      totalMapDist = atoi(optarg);
      break;
    case 'S':
      zygoteSelection = true;
      break;
    case 'm':
      selMult = atof(optarg);
      break;
    case 'e':
      selExp = atof(optarg);
      break;
    case 'r':
      seed = atoi(optarg);
      break;
    case 'o':
      outFile = optarg;
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
}

void getParameters( ostream & out,
                       int nbSimu,
                       int nbDiploids,
                       int nbGen,
                       int nbSitesPerChr,
                       int initNbTEsPerInd,
                       float probTransp0,
                       float k,
                       float probLoss,
                       int totalMapDist,
                       bool zygoteSelection,
                       float selMult,
                       float selExp,
                       int seed,
                       string outFile )
{
  out << "#nbSimu=" << nbSimu << endl;
  out << "#nbDiploids=" << nbDiploids << endl;
  out << "#nbGen=" << nbGen << endl;
  out << "#nbSitesPerChr=" << nbSitesPerChr << endl;
  out << "#initNbTEsPerInd=" << initNbTEsPerInd << endl;
  out << "#probTransp0=" << probTransp0 << endl;
  out << "#k=" << k << endl;
  out << "#probLoss=" << probLoss << endl;
  out << "#totalMapDist=" << totalMapDist << endl;
  out << "#zygoteSelection=" << boolalpha << zygoteSelection << noboolalpha << endl;
  out << "#selMult=" << selMult << endl;
  out << "#selExp=" << selExp << endl;
  out << "#seed=" << seed << endl;
  if( outFile != "" )
    out << "#output=" << outFile << endl;
}

void writeHeaderLine( ofstream & outStream )
{
  string sep = "\t";
  outStream << "simu" << sep << "gen"
            << sep << "nC" << sep << "meanC"
            << sep << "varC" << sep << "sdC"
            << sep << "minC" << sep << "q25C"
            << sep << "medC" << sep << "q75C"
            << sep << "maxC" << sep << "empty"
            << sep << "nL" << sep << "meanL"
            << sep << "varL" << sep << "sdL"
            << endl;
}

void getElapsedTime( ostream & out,
                     time_t startRawTime,
                     time_t endRawTime )
{
  time_t elapsedSec = difftime( endRawTime, startRawTime );
  tm * ptm;
  out << "#startTime: " << ctime( &startRawTime );
  out << "#endTime: " << ctime( &endRawTime );
  ptm = gmtime( &elapsedSec );
  out << "#elapsed time: "
      << setw(2) << setfill('0') << ptm->tm_hour << "h "
      << setw(2) << setfill('0') << ptm->tm_min << "m "
      << setw(2) << setfill('0') << ptm->tm_sec << "s"
      << endl;
}

int main( int argc, char* argv[] )
{
  int nbSimu = 1;
  int nbDiploids = 10;
  int nbGen = 10;
  int nbChrPerInd = 4;
  int nbSitesPerChr = 31;
  int initNbTEsPerInd = 10;
  float probLoss = 0.005;
  float probTransp0 = 0.01;
  int totalMapDist = 90;
  float k = 0.05;
  bool zygoteSelection = false;
  float selMult = 0.001;
  float selExp = 1.5;
  int seed = 1859;
  string outFile = "data.csv";
  int verbose = 0;
  gsl_rng * r;

  parse_args( argc, argv,
              nbSimu,
              nbDiploids,
              nbGen,
              nbSitesPerChr,
              initNbTEsPerInd,
              probTransp0,
              k,
              probLoss,
              totalMapDist,
              zygoteSelection,
              selMult,
              selExp,
              seed,
              outFile,
              verbose );

  time_t startRawTime;
  time( &startRawTime );
  printf ( "START: %s", ctime(&startRawTime) );

  if( verbose > 0 )
    getParameters( cout,
                   nbSimu,
                   nbDiploids,
                   nbGen,
                   nbSitesPerChr,
                   initNbTEsPerInd,
                   probTransp0,
                   k,
                   probLoss,
                   totalMapDist,
                   zygoteSelection,
                   selMult,
                   selExp,
                   seed,
                   outFile );

  // initialize outFile
  struct stat stFileInfo;
  int intStat;
  intStat = stat( outFile.c_str(), &stFileInfo );
  if( intStat == 0 )
    remove( outFile.c_str() );
  ofstream outStream;
  outStream.open( outFile.c_str(),
                  fstream::in | fstream::out | fstream::app );
  getParameters( outStream,
                 nbSimu,
                 nbDiploids,
                 nbGen,
                 nbSitesPerChr,
                 initNbTEsPerInd,
                 probTransp0,
                 k,
                 probLoss,
                 totalMapDist,
                 zygoteSelection,
                 selMult,
                 selExp,
                 seed,
                 "" );
  writeHeaderLine( outStream );

  // initialize the pseudo-random number generator
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc( T );
  gsl_rng_set( r, seed );

  // run the simulations
  for( int simuId=1; simuId<=nbSimu; ++simuId ){
    Simulation iSimu;
    iSimu.setSimulationIdentifier( simuId );
    iSimu.setNbGenerations( nbGen );
    iSimu.setNbDiploids( nbDiploids );
    iSimu.setNbChrPerIndividuals( nbChrPerInd );
    iSimu.setNbSitesPerChromosome( nbSitesPerChr );
    iSimu.setExpNbTEsPerIndividual( initNbTEsPerInd );
    iSimu.setTotalMapDist( totalMapDist );
    iSimu.setProbLoss( probLoss );
    iSimu.setProbTransp0( probTransp0 );
    iSimu.setK( k );
    iSimu.setZygoteSelection( zygoteSelection );
    iSimu.setSelMultiplicator( selMult );
    iSimu.setSelExponent( selExp );
    iSimu.setRng( r );
    iSimu.setOutFile( outFile );
    iSimu.setVerbose( verbose );
    iSimu.run();
  }

  gsl_rng_free( r );

  time_t endRawTime;
  time( &endRawTime );
  printf( "END: %s", ctime(&endRawTime) );

  getElapsedTime( outStream, startRawTime, endRawTime );
  outStream.close();
  if( verbose > 0 )
    getElapsedTime( cout, startRawTime, endRawTime );
}
