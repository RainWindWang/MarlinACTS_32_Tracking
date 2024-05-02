#include "ACTSKFTrackingProc.hxx"

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <LUXEEstimateTrackParamsFromSeed.hpp>
#include <LUXESeedfinder.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <Acts/TrackFitting/KalmanFitterError.hpp>

using namespace Acts::UnitLiterals;

#include "Helpers.hxx"
#include "MeasurementCalibrator.hxx"
#include "SeedSpacePoint.hxx"
#include "SourceLink.hxx"

// Track fitting definitions
// outlier to be looked at later on !!
using TrackFitterOptions =
    Acts::KalmanFitterOptions<ACTSTracking::MeasurementCalibrator,
                              ACTSTracking::VoidOutlierFinder>;
using TrackFitterResult =
      Acts::Result<Acts::KalmanFitterResult<ACTSTracking::SourceLink>>;

ACTSKFTrackingProc aACTSKFTrackingProc;

ACTSKFTrackingProc::ACTSKFTrackingProc()
    : ACTSProcBase("ACTSKFTrackingProc") {
  // modify processor description
  _description =
      "Fit tracks using the Kalman Filter algorithm";

  // Settings
  registerProcessorParameter(
      "InitialTrackError_RelP",
      "Track error estimate, momentum component (relative).",
      _initialTrackError_relP, 0.25);

  registerProcessorParameter("InitialTrackError_Phi",
                             "Track error estimate, phi (radians).",
                             _initialTrackError_phi, 1_degree);

  registerProcessorParameter("InitialTrackError_Lambda",
                             "Track error estimate, lambda (radians).",
                             _initialTrackError_lambda, 1_degree);

  registerProcessorParameter("InitialTrackError_Pos",
                             "Track error estimate, local position (mm).",
                             _initialTrackError_pos, 10_um);

  // Input collections - mc particles, tracker hits and the relationships
  // between them

  registerInputCollections(LCIO::TRACK, "TrackCollectionNames",
                           "Name of the Track input collections.",
                           _inputTrackCollections, {});

  registerOutputCollection(LCIO::TRACK, "TrackCollectionName",
                           "Name of track output collection.",
                           _outputTrackCollection, std::string("Tracks"));
}

void ACTSKFTrackingProc::init() {
  ACTSProcBase::init();

  // Reset counters
  _fitFails = 0;

}

void ACTSKFTrackingProc::processRunHeader(LCRunHeader *) {}

void ACTSKFTrackingProc::processEvent(LCEvent *evt) {
  //
  // Prepare the output
  // Make the output track collection
  LCCollectionVec *trackCollection = new LCCollectionVec(LCIO::TRACK);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  for (const std::string &collection : _inputTrackCollections) {
    LCCollection *trackCol = getCollection(collection, evt);
    if (trackCol == nullptr) continue;
  
    for (uint32_t idxTrk = 0;
         idxTrk < trackCol->getNumberOfElements(); idxTrk++) {
      EVENT::Track *trk = static_cast<EVENT::Track *>(
          trackCol->getElementAt(idxTrk));
      std::vector<EVENT::TrackerHit*> hits = trk->getTrackerHits();
      std::vector<std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *>> sortedHits;
      for (int i=0; i<hits.size(); i++){
        sortedHits.push_back(
          std::make_pair(geoIDMappingTool()->getGeometryID(hits[i]), hits[i]));
      }
      
  
    //
    // Prepare input hits in ACTS format

    // Loop over each hit collections and get a single vector with hits
    // from all of the subdetectors. Also include the Acts GeoId in
    // the vector. It will be important for the sort to speed up the
    // population of the final SourceLink multiset.
    // Sort by GeoID
    std::sort(
        sortedHits.begin(), sortedHits.end(),
        [](const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *> &hit0,
           const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *> &hit1)
            -> bool { return hit0.first < hit1.first; });

    // Turn the LCIO TrackerHit's into Acts objects
    // Assumes that the hits are sorted by the GeoID
    // Save first 3 hits as space points to later get initial parameters 
    // ACTSTracking::SourceLinkContainer sourceLinks;
    std::vector<ACTSTracking::SourceLink> sourceLinks;
    ACTSTracking::MeasurementContainer measurements;
    ACTSTracking::SeedSpacePointContainer spacePoints;

    sourceLinks.clear(); // irrelevant now but will be needed for multiple tracks
    sourceLinks.reserve(sortedHits.size());
    spacePoints.reserve(3);
    int i = 0;
    for (std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *> &hit :
         sortedHits) {
      // Convert to Acts hit
      const Acts::Surface *surface = trackingGeometry()->findSurface(hit.first);

      const double *lcioglobalpos = hit.second->getPosition();
      Acts::Vector3 globalPos = {lcioglobalpos[0], lcioglobalpos[1],
                                 lcioglobalpos[2]};
      Acts::Result<Acts::Vector2> lpResult =
          surface->globalToLocal(geometryContext(), globalPos, {0, 0, 0}, 0.5_um);
      //std::cout << surface << " " << trackingGeometry() << " " << globalPos[0] << " " << globalPos[2] << " " << hit.first << std::endl;
      if (!lpResult.ok())
        throw std::runtime_error(
            "Global to local transformation did not succeed.");

      Acts::Vector2 loc = lpResult.value();

      Acts::SymMatrix2 localCov = Acts::SymMatrix2::Zero();
      const EVENT::TrackerHitPlane *hitplane =
          dynamic_cast<const EVENT::TrackerHitPlane *>(hit.second);
      if (hitplane) {
        localCov(0, 0) = std::pow(hitplane->getdU() * Acts::UnitConstants::mm, 2);
        localCov(1, 1) = std::pow(hitplane->getdV() * Acts::UnitConstants::mm, 2);
      } else {
        throw std::runtime_error("Currently only support TrackerHitPlane.");
      }

      ACTSTracking::SourceLink sourceLink(surface->geometryId(),
                                          measurements.size(), hit.second);
      ACTSTracking::Measurement meas = Acts::makeMeasurement(
          sourceLink, loc, localCov, Acts::eBoundLoc0, Acts::eBoundLoc1);

      measurements.push_back(meas);
      // sourceLinks.emplace_hint(sourceLinks.end(), sourceLink);
      sourceLinks.push_back(sourceLink);

      //
      // Seed selection and conversion to useful coordinates
      if (i < 3) {
        Acts::RotationMatrix3 rotLocalToGlobal =
            surface->referenceFrame(geometryContext(), globalPos, {0, 0, 0});

        // Convert to a seed space point
        // the space point requires only the variance of the transverse and
        // longitudinal position. reduce computations by transforming the
        // covariance directly from local to rho/z.
        //
        // compute Jacobian from global coordinates to rho/z
        //
        //         rho = sqrt(x² + y²)
        // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
        //             = 2 * {x,y} / r
        //       dz/dz = 1 (duuh!)
        //
        double x = globalPos[Acts::ePos0];
        double y = globalPos[Acts::ePos1];
        double scale = 2 / std::hypot(x, y);
        Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
        jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
        jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
        jacXyzToRhoZ(1, Acts::ePos2) = 1;
        // compute Jacobian from local coordinates to rho/z
        Acts::ActsMatrix<2, 2> jac =
            jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
        // compute rho/z variance
        Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

        // Save spacepoint
        spacePoints.push_back(
            ACTSTracking::SeedSpacePoint(globalPos, var[0], var[1], sourceLink));
      }

      i++;
    }
    
    // save pointers of space points to spacePointPtrs vector
    std::vector<const ACTSTracking::SeedSpacePoint *> spacePointPtrs(
      spacePoints.size(), nullptr);
    std::transform(spacePoints.begin(), spacePoints.end(), spacePointPtrs.begin(),
                   [](const ACTSTracking::SeedSpacePoint &sp) { return &sp; });

    streamlog_out(DEBUG0) << "Created " << spacePoints.size() << " space points"
                          << std::endl;

    //
    // Caches
    Acts::MagneticFieldContext magFieldContext = Acts::MagneticFieldContext();
    Acts::MagneticFieldProvider::Cache magCache =
        magneticField()->makeCache(magFieldContext);

    //
    // Initialize track finder
    using Updater = Acts::GainMatrixUpdater;
    using Smoother = Acts::GainMatrixSmoother;
    using Stepper = Acts::EigenStepper<>;
    using Navigator = Acts::Navigator;
    using Propagator = Acts::Propagator<Stepper, Navigator>;
    using KF = Acts::KalmanFitter<Propagator, Updater, Smoother>;

    // Configurations
    Navigator::Config navigatorCfg{trackingGeometry()};
    navigatorCfg.resolvePassive = false;
    navigatorCfg.resolveMaterial = true;
    navigatorCfg.resolveSensitive = true;

    // construct all components for the fitter
    Stepper stepper(magneticField());
    Navigator navigator(navigatorCfg);
    Propagator propagator(std::move(stepper), std::move(navigator));
    KF trackFitter(std::move(propagator));
    
    Acts::PropagatorPlainOptions pOptions;
    pOptions.maxSteps = 10000;
    if (_propagateBackward) {
      pOptions.direction = Acts::backward;
    }

    // Construct a perigee surface as the target surface
    std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
        Acts::Surface::makeShared<Acts::PerigeeSurface>(
            Acts::Vector3{0., 0., 0.});

    std::unique_ptr<const Acts::Logger>
    logger=Acts::getDefaultLogger("TrackFitting",
    Acts::Logging::Level::INFO);//VERBOSE);

    // -- Outlier finder to be updated (currently accepts everything)
    TrackFitterOptions kfOptions = TrackFitterOptions(
        geometryContext(), magneticFieldContext(), calibrationContext(),
        ACTSTracking::MeasurementCalibrator(std::move(measurements)),
        ACTSTracking::VoidOutlierFinder{},
        Acts::LoggerWrapper{*logger}, pOptions,
        &(*perigeeSurface));

    const ACTSTracking::SeedSpacePoint *bottomSP = spacePointPtrs.front();
    const ACTSTracking::SourceLink &sourceLink = bottomSP->sourceLink();
    const Acts::GeometryIdentifier &geoId = sourceLink.geometryId();

    const Acts::Surface *surface = trackingGeometry()->findSurface(geoId);
    if (surface == nullptr) {
      std::cout << "surface with geoID " << geoId
                << " is not found in the tracking gemetry";
      // continue;
    }

    // Get the magnetic field at the bottom space point
    const Acts::Vector3 seedPos(bottomSP->x(), bottomSP->y(), bottomSP->z());
    Acts::Result<Acts::Vector3> seedField =
        magneticField()->getField(seedPos, magCache);
    if (!seedField.ok()) {
      throw std::runtime_error("Field lookup error: " +
                                seedField.error().value());
    }

    std::optional<Acts::BoundVector> optParams =
        Acts::estimateTrackParamsFromSeed(geometryContext(),
                                          spacePointPtrs.begin(), spacePointPtrs.end(),
                                          *surface, *seedField, 0.1_T);
    if (!optParams.has_value()) {
      std::cout << "Failed estimation of track parameters for seed."
                << std::endl;
      // continue;
    }

    const Acts::BoundVector &params = optParams.value();

    float charge = std::copysign(1, params[Acts::eBoundQOverP]);
    float p = std::abs(1 / params[Acts::eBoundQOverP]);

    // build the track covariance matrix using the smearing sigmas
    Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
        std::pow(_initialTrackError_pos, 2);
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
        std::pow(_initialTrackError_pos, 2);
    cov(Acts::eBoundTime, Acts::eBoundTime) =
        std::pow(_initialTrackError_time, 2);
    cov(Acts::eBoundPhi, Acts::eBoundPhi) =
        std::pow(_initialTrackError_phi, 2);
    cov(Acts::eBoundTheta, Acts::eBoundTheta) =
        std::pow(_initialTrackError_lambda, 2);
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
        std::pow(_initialTrackError_relP * p / (p * p), 2);

    Acts::BoundTrackParameters paramseed(surface->getSharedPtr(), params,
                                          charge, cov);

    streamlog_out(DEBUG) << "Seed Paramemeters" << std::endl
                          << paramseed << std::endl;

    //
    // Fit the track
    TrackFitterResult result = trackFitter.fit(sourceLinks, paramseed, kfOptions);

    if (result.ok()) {
        // Get the fit output object
        const Acts::KalmanFitterResult<ACTSTracking::SourceLink>
          &fitOutput = result.value();
        // The track entry indices container. One element here.
        std::vector<size_t> trackTips;
        trackTips.reserve(1);
        trackTips.emplace_back(fitOutput.lastMeasurementIndex);
        // The fitted parameters container. One element (at most) here.
        if (fitOutput.fittedParameters) {
          const Acts::BoundTrackParameters &params =
            fitOutput.fittedParameters.value();
          
          //
          // Helpful debug output
          Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
              Acts::MultiTrajectoryHelpers::trajectoryState(
                  fitOutput.fittedStates, trackTips[0]);
          streamlog_out(DEBUG) << "Trajectory Summary" << std::endl;
          streamlog_out(DEBUG)
              << "\tchi2Sum       " << trajState.chi2Sum << std::endl;
          streamlog_out(DEBUG)
              << "\tNDF           " << trajState.NDF << std::endl;
          streamlog_out(DEBUG)
              << "\tnHoles        " << trajState.nHoles << std::endl;
          streamlog_out(DEBUG)
              << "\tnMeasurements " << trajState.nMeasurements << std::endl;
          streamlog_out(DEBUG)
              << "\tnOutliers     " << trajState.nOutliers << std::endl;
          streamlog_out(DEBUG)
              << "\tnStates       " << trajState.nStates << std::endl;

          streamlog_out(DEBUG) << "Fitted Paramemeters" << std::endl
                        << params.parameters().transpose() << std::endl;
          
          // Make track object
          EVENT::Track *track = ACTSTracking::ACTS2Marlin_track(
              fitOutput, trackTips[0], magneticField(), magCache);

          // Save results
          trackCollection->addElement(track);

        } else {
          streamlog_out(WARNING) << "No fitted track parameters";
        }

      } else {
        streamlog_out(WARNING)
            << "Track fit error: " << result.error() << std::endl;
        _fitFails++;
      }
    }
  }
  // Save the output track collection
  evt->addCollection(trackCollection, _outputTrackCollection);
}

void ACTSKFTrackingProc::check(LCEvent *) {
  // nothing to check here - could be used to fill checkplots in reconstruction
  // processor
}

void ACTSKFTrackingProc::end() {}

LCCollection *ACTSKFTrackingProc::getCollection(
    const std::string &collectionName, LCEvent *evt) {
  try {
    return evt->getCollection(collectionName);
  } catch (DataNotAvailableException &e) {
    streamlog_out(DEBUG5) << "- cannot get collection. Collection "
                          << collectionName << " is unavailable" << std::endl;
    return nullptr;
  }
}
