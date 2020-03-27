/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2017 Albert Mink, Christopher McHardy
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#include "olb3D.h"
#include "olb3D.hh"   // include full template code

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;

typedef double T;

#define MINK
#define DESCRIPTOR descriptors::D3Q7<VELOCITY>

//#define DESCRIPTOR descriptors::D3Q27DescriptorRTLBM

SuperGeometry3D<T> prepareGeometryCube( T conversionFactorLength, T boundaryShift )
{
  OstreamManager clout( std::cout, "prepareGeometryCube" );
  clout << "Prepare cubeGeometry ..." << std::endl;

  Vector<T,3> origin = {0+boundaryShift,-0.5,-0.5};
  Vector<T,3> extent = {1-boundaryShift,1,1};
  IndicatorCuboid3D<T> span(extent, origin);
  origin = {0+boundaryShift,0,0};
  IndicatorCuboid3D<T> emittor(2*conversionFactorLength, 0.2, 0.2, origin);
  // infinite beam light source
  //IndicatorCuboid3D<T> emittor(2*conversionFactorLength, 1.0-2*conversionFactorLength, 1.0-2*conversionFactorLength, origin);

  CuboidGeometry3D<T>* cuboidGeometry = new CuboidGeometry3D<T>( span, conversionFactorLength, 2*singleton::mpi().getSize() );
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>( *cuboidGeometry );

  SuperGeometry3D<T> superGeometry( *cuboidGeometry, *loadBalancer, 2 );
  // set material number, 0 outside, 1 domain, 2 outflow, 3 inflow
  superGeometry.rename( 0, 2, span );
  superGeometry.rename( 2, 1, 1, 1, 1 );
  superGeometry.rename( 2, 3, emittor );
  // clean voxels with material number 0, inside and outside geometry
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare cubeGeometry ... OK" << std::endl;
  return superGeometry;
}

void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     T latticeRelaxationFrequency,
                     Dynamics<T,DESCRIPTOR>& bulkDynamics,
//                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bcDirected,
                     SuperGeometry3D<T>& superGeometry )
{
  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>(0) );
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );
  sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

  bcDirected.addTemperatureBoundary( superGeometry, 2, latticeRelaxationFrequency );
  bcDirected.addTemperatureBoundary( superGeometry, 3, latticeRelaxationFrequency );

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setBoundaryValues( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperGeometry3D<T>& superGeometry)
{
  // since adv-diffusion model is used, the velocity it set to 0
  AnalyticalConst3D<T,T> u0( 0.0, 0.0, 0.0 );        // 3D -> 3D
  AnalyticalConst3D<T,T> rho0( 0.0 );                // 3D -> 1D
  AnalyticalConst3D<T,T> rho1( 1.0 );                // 3D -> 1D

  // initialize media with density from analytical solution
  // at iT=0 the error is given by the maschinen genauigkeit
  sLattice.iniEquilibrium( superGeometry, 1, rho0, u0 );
  sLattice.iniEquilibrium( superGeometry, 2, rho0, u0 );
  sLattice.iniEquilibrium( superGeometry, 3, rho1, u0 );
  sLattice.defineRho( superGeometry, 2, rho0 );
  sLattice.defineRho( superGeometry, 3, rho1 );
  sLattice.initialize();
}

int main( int argc, char *argv[] )
{
  olbInit( &argc, &argv );
  string fName("app.xml");
  XMLreader config(fName);

  int RESOLUTION;
  T LATTICERELAXATIONTIME, ANISOTROPYFACTOR;
  config["Application"]["Discretization"]["Resolution"].read<int>(RESOLUTION);
  config["Application"]["Discretization"]["LatticeRelaxationTime"].read(LATTICERELAXATIONTIME);
  config["Application"]["PhysParameter"]["AnisotropyFactor"].read(ANISOTROPYFACTOR);

  // passed parameters might override resolution from XML file
//  std::string caseName{"case"};
//  if( argc < 3) {
//    caseName.append(argv[1]);
//  } else if ( argc < 4 ) {
//    caseName.append(argv[1]);
//    RESOLUTION = std::atoi(argv[2]);
//  }
  std::string caseName{"case1"};
#ifdef MINK
  singleton::directories().setOutputDir( "./"+std::to_string(RESOLUTION)+caseName+"_mink/" );
#else
  singleton::directories().setOutputDir( "./"+std::to_string(RESOLUTION)+caseName+"_mcHardy/" );
#endif
  OstreamManager clout( std::cout, "main" );
  clout << caseName << std::endl;


  T ABSORPTION, SCATTERING;
  std::stringstream xmlAbsorption(config["Application"]["PhysParameter"][caseName].getAttribute( "absorption"));
  xmlAbsorption >> ABSORPTION;
  std::stringstream xmlScattering(config["Application"]["PhysParameter"][caseName].getAttribute( "scattering"));
  xmlScattering >> SCATTERING;

  //LATTICERELAXATIONTIME = 1./(RESOLUTION*DESCRIPTOR<T>::invCs2*(1./(3*(ABSORPTION+SCATTERING))));
  LATTICERELAXATIONTIME = 1;// 1./(RESOLUTION*DESCRIPTOR<T>::invCs2*(1./(3*(ABSORPTION+SCATTERING)))+0.5);
  RadiativeUnitConverter<T,DESCRIPTOR> const converter( RESOLUTION, LATTICERELAXATIONTIME, ABSORPTION, SCATTERING, ANISOTROPYFACTOR );
#ifdef MINK
  clout << "working with diffuse approximation" << std::endl;
  T latticeSink = converter.getLatticeAbsorption()/converter.getLatticeDiffusion()/8.;

  PoissonDynamics<T,DESCRIPTOR> advDynamics( converter.getLatticeRelaxationTime(), instances::getBulkMomenta<T,DESCRIPTOR>(), latticeSink );
  // p1 motivated collision operator
  //P1Dynamics<T,DESCRIPTOR> advDynamics( LATTICERELAXATIONTIME, instances::getBulkMomenta<T,DESCRIPTOR>(), converter.getLatticeAbsorption(), converter.getLatticeScattering() );
#else
  clout << "working with direct discretization" << std::endl;
  int const q = descriptors::q<DESCRIPTOR>() - 1;
  std::array<std::array<double,q>, q> phi;
  double solution [q*(q +1)/2];
  computeAnisotropyMatrix<DESCRIPTOR>(1e-4, ANISOTROPYFACTOR, solution, phi);

  std::array<std::array<double,q+1>, q+1> anisoMatrix = {};
  for ( int m = 0; m < q; m++ ) {
    for ( int n = 0; n < q; n++ ) {
      anisoMatrix[m+1][n+1] = phi[m][n];
    }
  }
  RTLBMdynamicsMcHardyRK<T,DESCRIPTOR> advDynamics(
    instances::getBulkMomenta<T,DESCRIPTOR>(),
    converter.getLatticeAbsorption(),
    converter.getLatticeScattering(),
    anisoMatrix
    );
#endif
  converter.print();
  converter.write();

  SuperGeometry3D<T> superGeometry( prepareGeometryCube(converter.getConversionFactorLength(), 0. ) );
  //SuperGeometry3D<T> superGeometry( prepareGeometryCube(converter.getConversionFactorLength(), 1./converter.getExtinction() ) );
  SuperLattice3D<T,DESCRIPTOR> sLattice(superGeometry);
//  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryConditionDirected( sLattice );

  // addTemperatureBoundary results in Dirichlet Boundary for all directions, sum_i f_i = temperature
  // dynamics are AdvectionDiffusionBoundariesDynamics
//  createAdvectionDiffusionBoundaryCondition3D<T,DESCRIPTOR, AdvectionDiffusionBGKdynamics<T,DESCRIPTOR> >( sBoundaryCondition );

  // execute dynamics in src/dynamics/RtlbmBoundaryDynamics.h
  // inlet directed and fixed value (flux)
  createRtlbmDirectedBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryConditionDirected );
  // inlet diffuse and fixed value (flux)
  //createRtlbmDiffuseConstBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryConditionDirected );

  prepareLattice( sLattice,
                  converter.getLatticeRelaxationFrequency(),
                  advDynamics,
 //                 sBoundaryCondition,
                  sBoundaryConditionDirected,
                  superGeometry );
  int arbitrarySimulationTime = 13;
  Timer<double> timer( converter.getLatticeTime( arbitrarySimulationTime ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  setBoundaryValues(sLattice, converter, superGeometry);

  SuperVTMwriter3D<T> vtmWriter("data_app");
  SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice, superGeometry);
  vtmWriter.write( geometry );
  vtmWriter.createMasterFile();
  SuperLatticeDensity3D<T,DESCRIPTOR> density( sLattice );
  SuperLatticeFlux3D<T,DESCRIPTOR> flux( sLattice );
  vtmWriter.addFunctor(density);
  vtmWriter.addFunctor(flux);
  vtmWriter.addFunctor(geometry);

  util::ValueTracer<T> converge( 160, 1e-8);
  clout << "convergence criteria averages over: " << 160 << " iterations" << std::endl;

  int iT = 0;
  while( !(converge.hasConverged()) && iT<4*RESOLUTION )  {
    ++iT;
    sLattice.collideAndStream();
    converge.takeValue( sLattice.getStatistics().getAverageRho(), true );

    #ifdef MINK
    if ( iT % (2*RESOLUTION) == 0 ) {
    #else
    if ( iT % (RESOLUTION/5) == 0 ) {
    #endif
      sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
      timer.print( iT );
      timer.printStep();
      vtmWriter.write(iT);
    }
  }

  clout << "Simulation converged after " << iT << " iterations" << endl;
  AnalyticalFfromSuperF3D<T> analytDen( density, true, 1 );
  std::string str = std::to_string(singleton::directories().getLogOutDir()
                                   +std::to_string(RESOLUTION)+caseName+".csv");
  AnalyticalFfromSuperF3D<T> analytFlu( flux, true, 1 );
  clout << str << std::endl;
  FILE * pFile;
  const char * fileName = str.data();
  pFile  = fopen(fileName, "w");
  //fprintf(pFile, "%i\n", iT);
  //fprintf(pFile, "%s, %s, %s, %s\n", "position x", "0.0", "0.25", "0.375", "flux-x", "flux-y", "flux-z");
  for ( int nZ = 0; nZ <= 100; ++nZ ) {
    double line1[3] = {1.0*double(nZ)/100, 0, 0};
    double line2[3] = {1.0*double(nZ)/100, 0.25, 0};
    double line3[3] = {1.0*double(nZ)/100, 0.375, 0};
    double light1[1] = {0};
    double fluxx1[3] = {0.,0.,0.};
    double light2[1] = {0};
    double light3[1] = {0};
    analytDen( light1, line1 );
    analytFlu( fluxx1, line1 );
    analytDen( light2, line2 );
    analytDen( light3, line3 );
    if (singleton::mpi().getRank()==0) {
      printf("%4.3f, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", line1[0], light1[0], light2[0], light3[0],
        fluxx1[0], fluxx1[1], fluxx1[2]);
      fprintf(pFile, "%4.3f, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e, %8.6e\n", line1[0], light1[0], light2[0], light3[0],
        fluxx1[0], fluxx1[1], fluxx1[2]);
    }
  }
  fclose(pFile);

  timer.stop();
  timer.printSummary();

  return 0;
}
