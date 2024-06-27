#pragma once

#include <EVENT/LCEvent.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include "SourceLink.hxx"

namespace ACTSTracking {

using TrackResult = Acts::TrackContainer<Acts::VectorTrackContainer,
                                         Acts::VectorMultiTrajectory,
                                         std::shared_ptr>::TrackProxy;

/// @brief void outlier finder
// an outlier finder that accepts any hit (does not reject anything)
struct VoidOutlierFinder {
/// @brief Public call mimicking an outlier finder
///
/// @tparam track_state_t Type of the track state
///
/// @param trackState The trackState to investigate
///
/// @return Whether it's outlier or not
template <typename track_state_t>
constexpr bool operator()(const track_state_t& trackState) const {
(void)trackState;
return false;
}
};

//! Get path to a resource file
/**
 * Get absolute file of a file `inpath` by looking in the following places:
 *  - `inpath` to the current working directory
 *  - `ACTSTRACKING_SOURCEDIR/inpath`
 *  - `ACTSTRACKING_DATADIR/inpath`
 *
 * If the files is not found at any location, then `inpath` is returned.
 * If `path` starts with a /, then it is returned directly.
 *
 * \parm inpath File to find.
 *
 * \return Absolute path to file.
 */
std::string findFile(const std::string& inpath);


/*
// this was used in old ACTS2Marlin_track
inline auto getFitParams(const Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>& fitOutput, std::size_t trackTip) {
    return fitOutput.fittedParameters.at(trackTip);
}

inline auto getFitParams(const Acts::KalmanFitterResult<ACTSTracking::SourceLink>& fitOutput, std::size_t) {
    return fitOutput.fittedParameters.value();
}
*/
//! Convert ACTS CKF result to LCIO track class -------------------------------------------------------------------------
/**
 * Converted propertie are:
 *  - goodness of fit (chi2, ndf)
 *  - associated hits
 *  - track states at IP
 *
 * \param fitOutput CKF fit result
 * \param trackTip index of track to convert inside fit result
 * \param magneticField magnetic field at different locations in the detector
 * \param magCache cache to help with magnetic field lookup
 *
 * \return Track with equivalent parameters of the ACTS track
 */
 /*
 template<typename KFResultT>
EVENT::Track* ACTS2Marlin_track(
    const KFResultT&
        fitOutput,
    std::size_t trackTip,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache) {
          IMPL::TrackImpl* track = new IMPL::TrackImpl;
  Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates,
                                                    trackTip);
  // Fit state
  track->setChi2(trajState.chi2Sum);
  track->setNdf(trajState.NDF);

  // Track state at IP -------------------------------------------------------------------------
  static const Acts::Vector3 zeroPos(0, 0, 0);
  Acts::Result<Acts::Vector3> fieldRes =
      magneticField->getField(zeroPos, magCache);
  if (!fieldRes.ok()) {
    throw std::runtime_error("Field lookup error: " + fieldRes.error().value());
  }
  Acts::Vector3 field = *fieldRes;

  const Acts::BoundTrackParameters& params = getFitParams(fitOutput, trackTip);

  EVENT::TrackState* trackStateAtIP = ACTSTracking::ACTS2Marlin_trackState(
      EVENT::TrackState::AtIP, params, field[2] / Acts::UnitConstants::T);
  track->trackStates().push_back(trackStateAtIP);

  //
  // Hits and associated track states
  EVENT::TrackerHitVec hitsOnTrack;
  EVENT::TrackStateVec statesOnTrack;
  fitOutput.fittedStates.visitBackwards(
      trackTip, [&](const Acts::MultiTrajectory<
                    ACTSTracking::SourceLink>::ConstTrackStateProxy& state) {
        // No measurement at this state
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        // register all particles that generated this hit
        EVENT::TrackerHit* myHit = state.uncalibrated().lciohit();
        hitsOnTrack.push_back(myHit);

        // Save track state information
        const Acts::Vector3 hitPos(myHit->getPosition()[0],
                                   myHit->getPosition()[1],
                                   myHit->getPosition()[2]);
        Acts::Result<Acts::Vector3> fieldRes =
            magneticField->getField(hitPos, magCache);
        if (!fieldRes.ok()) {
          throw std::runtime_error("Field lookup error: " +
                                   fieldRes.error().value());
        }
        Acts::Vector3 field = *fieldRes;

        EVENT::TrackState* trackState = ACTSTracking::ACTS2Marlin_trackState(
            EVENT::TrackState::AtOther, state.smoothed(),
            state.smoothedCovariance(), field[2] / Acts::UnitConstants::T);
        statesOnTrack.push_back(trackState);

        return true;
      });

  // Reverse hits and states, above creates them backwards
  std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
  std::reverse(statesOnTrack.begin(), statesOnTrack.end());

  // Save hits
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(
      lcio::LCTrackerCellID::encoding_string());
  EVENT::IntVec& subdetectorHitNumbers = track->subdetectorHitNumbers();
  for (EVENT::TrackerHit* hit : hitsOnTrack) {
    track->addHit(hit);

    uint32_t sysid = decoder(hit)["system"];
    if (subdetectorHitNumbers.size() <= sysid) {
      subdetectorHitNumbers.resize(sysid + 1, 0);
    }
    subdetectorHitNumbers[sysid]++;
  }

  // Save the track states at hits
  if (statesOnTrack.size() > 0) {
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.back())
        ->setLocation(EVENT::TrackState::AtLastHit);
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.front())
        ->setLocation(EVENT::TrackState::AtFirstHit);
  }

  EVENT::TrackStateVec& myTrackStates = track->trackStates();
  myTrackStates.insert(myTrackStates.end(), statesOnTrack.begin(),
                       statesOnTrack.end());

  return track;
}*/

// EVENT::Track* ACTS2Marlin_track_KF( 
//     const Acts::KalmanFitterResult<ACTSTracking::SourceLink>&
//         fitOutput,
//     std::size_t trackTip,
//     std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
//     Acts::MagneticFieldProvider::Cache& magCache);

//! Convert ACTS KF result to LCIO track class
/**
 * Converted propertie are:
 *  - goodness of fit (chi2, ndf)
 *  - associated hits
 *  - track states at IP
 *
 * \param fitOutput KF fit result
 * \param magneticField magnetic field at different locations in the detector
 * \param magCache cache to help with magnetic field lookup
 *
 * \return Track with equivalent parameters of the ACTS track
 */

//! Convert ACTS result to LCIO track class
EVENT::Track* ACTS2Marlin_track(
    const TrackResult& fitter_res, // TrackResult instead of KFResult & CKFResult
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache);

//! Convert ACTS track state class to Marlin class -------------------------------------------------------------------------
/**
 * \param location Location where the track state is defined (ie: `AtIP`)
 * \param params ACTS track state parameters
 * \params Bz magnetic field at location of track state [Tesla]
 *
 * \return Track state with equivalent parameters of the ACTS track
 */
EVENT::TrackState* ACTS2Marlin_trackState(
    int location, const Acts::BoundTrackParameters& params, double Bz);

EVENT::TrackState* ACTS2Marlin_trackState(int location,
                                          const Acts::BoundVector& value,
                                          const Acts::BoundMatrix& cov,
                                          double Bz);

//! Get collection from `LCEvent` with silent fail
/**
 * \param evt event store
 * \param event collection name
 *
 * \return Collection, if found, `nullptr` otherwise
 */
EVENT::LCCollection* getCollection(EVENT::LCEvent* evt,
                                   const std::string& name);

Acts::ParticleHypothesis convertParticle(const EVENT::MCParticle* mcParticle);

}  // namespace ACTSTracking
