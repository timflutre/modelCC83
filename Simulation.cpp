/*
 * \file Simulation.cpp
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
using namespace std;

#include "Simulation.h"
#include "Population.h"

Simulation::Simulation( void )
{
  setSimulationIdentifier( 0 );
  setNbGenerations( 0 );
  setNbDiploids( 0 );
  setNbChrPerIndividuals( 0 );
  setNbSitesPerChromosome( 0 );
  setExpNbTEsPerIndividual( 0 );
  setTotalMapDist( 0 );
  setProbLoss( 0.0 );
  setProbTransp0( 0.0 );
  setK( 0.0 );
  setZygoteSelection( false );
  setSelMultiplicator( 0.0 );
  setSelExponent( 0.0 );
  setOutFile( "data.tsv" );
  setVerbose( 0 );
}

void Simulation::setSimulationIdentifier( int si )
{
  simuId = si;
}

void Simulation::setNbGenerations( int ng )
{
  nbGen = ng;
}

void Simulation::setNbDiploids( int nd )
{
  nbDiploids = nd;
}

void Simulation::setNbChrPerIndividuals( int cpi )
{
  nbChrPerInd = cpi;
}

void Simulation::setNbSitesPerChromosome( int spc )
{
  nbSitesPerChr = spc;
}

void Simulation::setExpNbTEsPerIndividual( int nti )
{
  expNbTEsPerInd = nti;
}

void Simulation::setTotalMapDist( int tmd )
{
  totalMapDist = tmd;
}

void Simulation::setProbLoss( float pl )
{
  probLoss = pl;
}

void Simulation::setProbTransp0( float pt )
{
  probTransp0 = pt;
}

void Simulation::setK( float param_k )
{
  k = param_k;
}

void Simulation::setZygoteSelection( bool zs )
{
  zygoteSelection = zs;
}

void Simulation::setSelMultiplicator( float sm )
{
  selMult = sm;
}

void Simulation::setSelExponent( float se )
{
  selExp = se;
}

void Simulation::setOutFile( string of )
{
  outFile = of;
}

void Simulation::setVerbose( int v )
{
  verbose = v;
}

void Simulation::setRng( gsl_rng * rng )
{
  r = rng;
}

int Simulation::getSimulationIdentifier( void )
{
  return( simuId );
}

int Simulation::getNbGenerations( void )
{
  return( nbGen );
}

int Simulation::getNbDiploids( void )
{
  return( nbDiploids );
}

int Simulation::getNbChrPerIndividual( void )
{
  return( nbChrPerInd );
}

int Simulation::getNbSitesPerChromosome( void )
{
  return( nbSitesPerChr );
}

int Simulation::getExpNbTEsPerIndividual( void )
{
  return( expNbTEsPerInd );
}

int Simulation::getTotalMapDist( void )
{
  return( totalMapDist );
}

float Simulation::getProbLoss( void )
{
  return( probLoss );
}

float Simulation::getProbTransp0( void )
{
  return( probTransp0 );
}

float Simulation::getK( void )
{
  return( k );
}

bool Simulation::getZygoteSelection( void )
{
  return( zygoteSelection );
}

float Simulation::getSelMultiplicator( void )
{
  return( selMult );
}

float Simulation::getSelExponent( void )
{
  return( selExp );
}

string Simulation::getOutFile( void )
{
  return( outFile );
}

int Simulation::getVerbose( void )
{
  return( verbose );
}

void Simulation::printSimGen( int g )
{
  cout << "simulation " << simuId
       << ": generation " << setw(4) << setfill('0') << g
       << "/" << nbGen << endl;;
}

void Simulation::run( void )
{
  Population pop;
  pop.setNbDiploids( getNbDiploids() );
  pop.setNbChrPerIndividual( getNbChrPerIndividual() );
  pop.setNbSitesPerChromosome( getNbSitesPerChromosome() );
  pop.setExpNbTEsPerIndividual( getExpNbTEsPerIndividual() );
  pop.setTotalMapDist( getTotalMapDist() );
  pop.setZygoteSelection( getZygoteSelection() );
  pop.setSelMultiplicator( getSelMultiplicator() );
  pop.setSelExponent( getSelExponent() );
  pop.setVerbose( getVerbose()-1 );
  pop.setRng( r );
  pop.initialize();
  pop.saveData( getSimulationIdentifier(),
                0, getOutFile() );

  for( int g=1; g<=nbGen; ++g ){
    if( getVerbose() > 0 ){
      printSimGen( g );
      pop.printDistribTEsPerInd();
    }
    if( pop.getSumNbTEs() > 0 ){
      pop.makeNewGeneration( getVerbose()-1 );
      pop.loss( probLoss );
      pop.transposition( probTransp0, k );
      pop.saveData( getSimulationIdentifier(),
                    g, getOutFile() );
    }
    else
      break;
  }
}
