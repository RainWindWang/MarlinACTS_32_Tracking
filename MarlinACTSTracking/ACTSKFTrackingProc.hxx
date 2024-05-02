#ifndef ACTSKFTrackingProc_h
#define ACTSKFTrackingProc_h 1

#include <EVENT/TrackerHit.h>

#include <UTIL/CellIDDecoder.h>

#include <Acts/Definitions/Units.hpp>

#include "ACTSProcBase.hxx"
#include "GeometryIdSelector.hxx"

/**
 * This code performs a true pattern recognition by looping over all MC
 * particles and adding all hits associated to them onto a prototrack. This is
 * then fitted and output.
 */
class ACTSKFTrackingProc : public ACTSProcBase {
 public:
  virtual marlin::Processor *newProcessor() {
    return new ACTSKFTrackingProc;
  }

  ACTSKFTrackingProc(const ACTSKFTrackingProc &) = delete;
  ACTSKFTrackingProc &operator=(const ACTSKFTrackingProc &) =
      delete;
  ACTSKFTrackingProc();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

 private:
  /** Call to get collections
   */
  LCCollection *getCollection(const std::string &, LCEvent *);

 protected:
  // Collection names for (in/out)put
  std::vector<std::string> _inputTrackCollections;
  std::string _outputTrackCollection;

  // Run settings
  bool _propagateBackward = false;

  // Track fit parameters
  double _initialTrackError_pos;
  double _initialTrackError_phi;
  double _initialTrackError_relP;
  double _initialTrackError_lambda;
  double _initialTrackError_time =
      100 * Acts::UnitConstants::ns;  // No Marlin default

  uint32_t _fitFails;
};

#endif
