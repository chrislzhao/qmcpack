//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_QMCDRIVERNEWTESTWRAPPER_H
#define QMCPLUSPLUS_QMCDRIVERNEWTESTWRAPPER_H
#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/DriverTraits.h"

namespace qmcplusplus
{
/** Unit testing an impure virtual base class
 *  requires a absolute minimal subtype
 */
class QMCDriverNewTestWrapper : public QMCDriverNew
{
public:
  QMCDriverNewTestWrapper(QMCDriverInput&& input,
                          MCPopulation& population,
                          TrialWaveFunction& psi,
                          QMCHamiltonian& h,
                          WaveFunctionPool& ppool,
                          Communicate* comm)
      : QMCDriverNew(std::move(input), population, psi, h, ppool, comm)
  {}

  QMCRunType getRunType() { return QMCRunType::DUMMY; }
  IndexType calc_default_local_walkers() { return 0; };
  bool run() { return false; }
};

} // namespace qmcplusplus
#endif
